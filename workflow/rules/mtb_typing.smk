rule mtb_lineage_id:
    input:
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
    output:
        tsv=OUT + "/mtb_typing/lineage_call/{sample}.tsv",
    conda:
        "../envs/fast_lineage_caller.yaml"
    log:
        OUT + "/log/mtb_lineage_id/{sample}.log",
    message:
        "Typing Mtb lineage for {wildcards.sample}"
    threads: 1
    resources:
        mem_gb=config["mem_gb"]["other"],
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
    log:
        OUT + "/log/mtb_typing/mtb_coll_contamination/{sample}.log",
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
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
    output:
        tsv=OUT + "/mtb_typing/contamination_check/rrs_rrl_contamination/{sample}.tsv",
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/mtb_typing/mtb_rrs_rrl_contamination/{sample}.log",
    shell:
        """
gatk CountVariants \
-V {input.vcf} \
1> {output.tsv} \
2> {log}
        """


rule mtb_snpeff_annotation:
    input:
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
        db_dir=OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff_ref",
        config=OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff.config",
        dummy=OUT + "/mtb_typing/prepared_reference_data/{sample}/build_snpeff_db.done",
    output:
        vcf=OUT + "/mtb_typing/snpeff/{sample}.vcf",
        stats=OUT + "/mtb_typing/snpeff_stats/{sample}.html",
    conda:
        "../envs/snpeff.yaml"
    # params:
    #     outdir = OUT,
    #     db_dir_base = "snpeff_ref"
    log:
        OUT + "/log/mtb_snpeff_annotation/{sample}.log",
    shell:
        """
WORKDIR=$(dirname {input.db_dir})
DB_NAME=$(basename {input.db_dir})
snpEff ann -c {input.config} -dataDir $WORKDIR -stats {output.stats} -noLog -o gatk -ud 0 $DB_NAME {input.vcf} > {output.vcf} 2>{log}
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
    params:
        base_annotations="CHROM,POS",
        extra_annotations=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_columns"
        ],
    log:
        OUT + "/log/mtb_annotate_ab_positions/{sample}.log",
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
        tsv=OUT + "/mtb_typing/annotated_variants/{sample}.tsv",
    conda:
        "../envs/gatk_picard.yaml"
    params:
        metadata=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_columns"
        ],
    log:
        OUT + "/log/mtb_annotated_vcf_to_table/{sample}.log",
    shell:
        """
FIELDS=$(python workflow/scripts/print_fields_VariantsToTable.py {params.metadata:q})
gatk VariantsToTable \
-V {input.vcf} \
-F CHROM \
-F POS \
-F TYPE \
-F REF \
-F ALT \
-F DP \
-GF AF \
$FIELDS \
-O {output.tsv} 2>&1>{log}
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
