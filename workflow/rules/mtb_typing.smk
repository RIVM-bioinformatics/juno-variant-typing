rule mtb_lineage_id:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.vcf",
    output:
        tsv=OUT + "/mtb_typing/lineage_call/{sample}.tsv",
    conda:
        "../envs/fast_lineage_caller.yaml"
    container:
        "docker://ghcr.io/boasvdp/fast_lineage_caller:1.0.0"
    log:
        OUT + "/log/mtb_lineage_id/{sample}.log",
    message:
        "Typing Mtb lineage for {wildcards.sample}"
    threads: config["threads"]["fast-lineage-caller"]
    resources:
        mem_gb=config["mem_gb"]["fast-lineage-caller"],
    shell:
        """
fast-lineage-caller \
--out {output.tsv} \
{input.vcf} \
2>&1>{log}
        """


rule mtb_coll_contamination:
    input:
        bam=OUT + "/mtb_typing/prepared_files/{sample}.bam",
        bai=OUT + "/mtb_typing/prepared_files/{sample}.bam.bai",
        reference=OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta",
        dummy=OUT + "/mtb_typing/prepared_files/{sample}_ref.dict",
        fai=OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta.fai",
        bed=lambda wildcards: SAMPLES[wildcards.sample]["single_copy_bed"],
    output:
        tsv=OUT + "/mtb_typing/contamination_check/coll_positions/{sample}.tsv",
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/mtb_typing/mtb_coll_contamination/{sample}.log",
    message:
        "Assessing minor variants in coll positions for {wildcards.sample}"
    threads: config["threads"]["gatk"]
    resources:
        mem_gb=config["mem_gb"]["gatk"],
    shell:
        """
gatk CollectAllelicCounts \
-I {input.bam} \
-R {input.reference} \
-L {input.bed} \
-O {output.tsv} \
2>&1>{log}
        """


rule mtb_rrs_rrl_contamination:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.vcf",
        tbi=OUT + "/mtb_typing/prepared_files/{sample}.vcf.idx",
        bed=lambda wildcards: SAMPLES[wildcards.sample]["count_mutations_bed"],
    output:
        filtered_vcf=temp(
            OUT + "/mtb_typing/contamination_check/rrs_rrl_contamination/{sample}.vcf"
        ),
        tsv=OUT + "/mtb_typing/contamination_check/rrs_rrl_contamination/{sample}.tsv",
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/mtb_typing/mtb_rrs_rrl_contamination/{sample}.log",
    message:
        "Counting rrs/rrl variants for {wildcards.sample}"
    threads: config["threads"]["gatk"]
    resources:
        mem_gb=config["mem_gb"]["gatk"],
    shell:
        """
gatk SelectVariants \
-V {input.vcf} \
-O {output.filtered_vcf} \
--exclude-filtered \
2>&1>{log}

gatk CountVariants \
-V {output.filtered_vcf} \
-L {input.bed} \
1> {output.tsv} \
2>> {log}
        """


rule mtb_snpeff_annotation:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.vcf",
        db_dir=OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff_ref",
        config=OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff.config",
        dummy=OUT + "/mtb_typing/prepared_reference_data/{sample}/build_snpeff_db.done",
    output:
        vcf=OUT + "/mtb_typing/snpeff/{sample}.vcf",
        stats=OUT + "/mtb_typing/snpeff_stats/{sample}.html",
    conda:
        "../envs/snpeff.yaml"
    container:
        "docker://staphb/snpeff:5.1"
    log:
        OUT + "/log/mtb_snpeff_annotation/{sample}.log",
    params:
        use_singularity=config["use_singularity"],
    message:
        "Running SnpEff for {wildcards.sample}"
    threads: config["threads"]["snpeff"]
    resources:
        mem_gb=config["mem_gb"]["snpeff"],
    shell:
        """
if [ {params.use_singularity} == True ]
then
    EXEC=snpeff
else
    EXEC=snpEff
fi
WORKDIR=$(dirname {input.db_dir})
DB_NAME=$(basename {input.db_dir})
$EXEC ann -c {input.config} -dataDir $WORKDIR -stats {output.stats} -noLog -o gatk -ud 0 $DB_NAME {input.vcf} > {output.vcf} 2>{log}
        """


rule mtb_annotated_vcf_to_table:
    input:
        vcf=OUT + "/mtb_typing/snpeff/{sample}.vcf",
    output:
        tsv=temp(OUT + "/mtb_typing/annotated_variants/raw/{sample}.tsv"),
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    params:
        metadata=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_columns"
        ],
        effect_column="EFF",
    log:
        OUT + "/log/mtb_annotated_vcf_to_table/{sample}.log",
    message:
        "Convert annotated variants to table for {wildcards.sample}"
    threads: config["threads"]["gatk"]
    resources:
        mem_gb=config["mem_gb"]["gatk"],
    shell:
        """
gatk VariantsToTable \
-V {input.vcf} \
--show-filtered \
-F CHROM \
-F POS \
-F TYPE \
-F REF \
-F ALT \
-F DP \
-F FILTER \
-GF AF \
-F {params.effect_column} \
-O {output.tsv} 2>&1>{log}
        """


rule mtb_annotate_ab_positions:
    input:
        tsv=OUT + "/mtb_typing/annotated_variants/raw/{sample}.tsv",
        reslist=lambda wildcards: SAMPLES[wildcards.sample]["resistance_variants_csv"],
    output:
        tsv=OUT + "/mtb_typing/annotated_variants/{sample}.tsv",
    params:
        merge_cols="CHROM,POS",
        keep_cols=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_columns"
        ],
    log:
        OUT + "/log/mtb_annotate_ab_positions/{sample}.log",
    message:
        "Annotating variants with AMR for {wildcards.sample}"
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
python workflow/scripts/postprocess_variant_table.py \
--input {input.tsv} \
--reference_data {input.reslist} \
--merge_cols {params.merge_cols} \
--keep_cols {params.keep_cols} \
--output {output} \
2>&1>{log}
        """


rule mtb_filter_res_table_positions:
    input:
        tsv=OUT + "/mtb_typing/annotated_variants/{sample}.tsv",
    output:
        tsv=OUT + "/mtb_typing/annotated_resistance_filtered/{sample}.tsv",
    params:
        ab_column=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_ab_column"
        ],
    log:
        OUT + "/log/mtb_filter_res_table_positions/{sample}.log",
    shell:
        """
python workflow/scripts/filter_res_table.py \
--input {input} \
--ab-column {params.ab_column} \
--output {output} 2>&1>{log}
        """


rule mtb_make_json:
    input:
        lineage=OUT + "/mtb_typing/lineage_call/{sample}.tsv",
        rrs_rrl_snp_counts=OUT
        + "/mtb_typing/contamination_check/rrs_rrl_contamination/{sample}.tsv",
        picard_collectwgsmetrics=INPUT + "/qc_mapping/CollectWgsMetrics/{sample}.txt",
    output:
        json=OUT + "/mtb_typing/seq_exp_json/{sample}.json",
    log:
        OUT + "/log/mtb_make_json/{sample}.log",
    shell:
        """
python workflow/scripts/create_tb_json.py \
--lineage-call {input.lineage} \
--rrs-rrl-snp-counts {input.rrs_rrl_snp_counts} \
--picard {input.picard_collectwgsmetrics} \
--output {output.json} 2>&1>{log}
        """


module consensus_workflow:
    config:
        config
    snakefile:
        "consensus.smk"


use rule mark_variants_by_proximity from consensus_workflow with:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.vcf",
    output:
        bed=OUT + "/mtb_typing/prepared_files/{sample}.no_prox.bed",
    log:
        OUT + "/log/mark_variants_by_proximity/{sample}.log",


use rule subset_fixed_snps_from_vcf from consensus_workflow with:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.vcf",
    output:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.fixed_snps.vcf",
    log:
        OUT + "/log/subset_fixed_snps_from_vcf/{sample}.log",


use rule subset_low_confidence_variants_from_vcf from consensus_workflow with:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.vcf",
    output:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.lowconf.vcf",
    log:
        OUT + "/log/subset_low_confidence_variants_from_vcf/{sample}.log",


use rule zip_and_index_sample_vcf_bcftools from consensus_workflow with:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.fixed_snps.vcf",
    output:
        vcf_gz=OUT + "/mtb_typing/prepared_files/{sample}.fixed_snps.vcf.gz",
        tbi=OUT + "/mtb_typing/prepared_files/{sample}.fixed_snps.vcf.gz.tbi",
    log:
        OUT + "/log/zip_and_index_sample_vcf_bcftools/{sample}.log",


use rule introduce_mutations_to_reference from consensus_workflow with:
    input:
        vcf_gz=OUT + "/mtb_typing/prepared_files/{sample}.fixed_snps.vcf.gz",
        tbi=OUT + "/mtb_typing/prepared_files/{sample}.fixed_snps.vcf.gz.tbi",
        reference=OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta",
    output:
        fasta=OUT + "/mtb_typing/consensus/raw/{sample}.fasta",
    log:
        OUT + "/log/introduce_mutations_to_reference/{sample}.log",


use rule mask_fasta_on_depth_from_bam from consensus_workflow with:
    input:
        bam=OUT + "/mtb_typing/prepared_files/{sample}.bam",
        bai=OUT + "/mtb_typing/prepared_files/{sample}.bam.bai",
        fasta=OUT + "/mtb_typing/consensus/raw/{sample}.fasta",
    output:
        fasta=OUT + "/mtb_typing/consensus/depth_masked/{sample}.fasta",
    log:
        OUT + "/log/mask_fasta_on_depth_from_bam/{sample}.log",


use rule mask_fasta_based_on_bed_or_vcf from consensus_workflow as mask_fasta_on_low_confidence_variants with:
    input:
        features=OUT + "/mtb_typing/prepared_files/{sample}.lowconf.vcf",
        fasta=OUT + "/mtb_typing/consensus/depth_masked/{sample}.fasta",
    output:
        fasta=OUT + "/mtb_typing/consensus/depth_masked_low_conf_masked/{sample}.fasta",
    log:
        OUT + "/log/mask_fasta_on_low_confidence_variants/{sample}.log",


use rule mask_fasta_based_on_bed_or_vcf from consensus_workflow as mask_fasta_on_proximity_variants with:
    input:
        features=OUT + "/mtb_typing/prepared_files/{sample}.no_prox.bed",
        fasta=OUT + "/mtb_typing/consensus/depth_masked_low_conf_masked/{sample}.fasta",
    output:
        fasta=OUT
        + "/mtb_typing/consensus/depth_masked_low_conf_masked_proxmask/{sample}.fasta",
    log:
        OUT + "/log/mask_fasta_on_proximity_variants/{sample}.log",


if Path(INPUT + "/variants_raw/mask.bed").is_file():

    use rule mask_fasta_based_on_bed_or_vcf from consensus_workflow as mask_fasta_by_custom_bed with:
        input:
            features=INPUT + "/variants_raw/mask.bed",
            fasta=OUT
            + "/mtb_typing/consensus/depth_masked_low_conf_masked_proxmask/{sample}.fasta",
        output:
            fasta=OUT
            + "/mtb_typing/consensus/depth_masked_low_conf_masked_proxmask_bedmasked/{sample}.fasta",
        log:
            OUT + "/log/mask_fasta_on_custom_bed/{sample}.log",

else:

    rule:
        input:
            fasta=OUT
            + "/mtb_typing/consensus/depth_masked_low_conf_masked_proxmask/{sample}.fasta",
        output:
            fasta=OUT
            + "/mtb_typing/consensus/depth_masked_low_conf_masked_proxmask_bedmasked/{sample}.fasta",
        shell:
            """
cp {input} {output}
            """


use rule replace_fasta_header from consensus_workflow with:
    input:
        fasta=OUT
        + "/mtb_typing/consensus/depth_masked_low_conf_masked_proxmask_bedmasked/{sample}.fasta",
    output:
        fasta=OUT + "/mtb_typing/consensus/{sample}.fasta",
    log:
        OUT + "/log/replace_fasta_header/{sample}.log",
