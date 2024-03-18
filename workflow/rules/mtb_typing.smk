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


rule mtb_annotate_ab_positions:
    input:
        vcf=OUT + "/mtb_typing/snpeff/{sample}.vcf",
        compressed_table=OUT
        + "/mtb_typing/prepared_reference_data/{sample}/ab_table.tab.gz",
        header=OUT + "/mtb_typing/prepared_reference_data/{sample}/ab_table.header",
    output:
        vcf=OUT + "/mtb_typing/annotated_vcf/{sample}.vcf",
    conda:
        "../envs/bcftools.yaml"
    container:
        "docker://staphb/bcftools:1.19"
    params:
        base_annotations="CHROM,POS",
        extra_annotations=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_columns"
        ],
    log:
        OUT + "/log/mtb_annotate_ab_positions/{sample}.log",
    message:
        "Annotating variants with AMR for {wildcards.sample}"
    threads: config["threads"]["bcftools"]
    resources:
        mem_gb=config["mem_gb"]["bcftools"],
    shell:
        """
bcftools annotate \
-a {input.compressed_table} \
-h {input.header} \
-c {params.base_annotations},{params.extra_annotations:q} \
{input.vcf} \
1> {output.vcf} \
2> {log}
        """


rule mtb_annotated_vcf_to_table:
    input:
        vcf=OUT + "/mtb_typing/annotated_vcf/{sample}.vcf",
    output:
        tsv=temp(OUT + "/mtb_typing/annotated_variants/raw_{sample}.tsv"),
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
FIELDS=$(python workflow/scripts/print_fields_VariantsToTable.py {params.metadata:q},{params.effect_column:q})
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
$FIELDS \
-O {output.tsv} 2>&1>{log}
        """


rule postprocess_variant_table:
    input:
        tsv=OUT + "/mtb_typing/annotated_variants/raw_{sample}.tsv",
    output:
        tsv=OUT + "/mtb_typing/annotated_variants/{sample}.tsv",
    log:
        OUT + "/log/postprocess_variant_table/{sample}.log",
    shell:
        """
python workflow/scripts/postprocess_variant_table.py \
--input {input} \
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
