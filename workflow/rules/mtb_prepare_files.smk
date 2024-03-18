rule copy_sample_bam:
    input:
        bam=lambda wildcards: SAMPLES[wildcards.sample]["bam"],
    output:
        bam=temp(OUT + "/mtb_typing/prepared_files/{sample}.bam"),
    log:
        OUT + "/log/copy_sample_bam/{sample}.log",
    shell:
        """
cp {input.bam} {output.bam} 2>&1>{log}
        """


rule index_sample_bam:
    input:
        bam=OUT + "/mtb_typing/prepared_files/{sample}.bam",
    output:
        bai=temp(OUT + "/mtb_typing/prepared_files/{sample}.bam.bai"),
    container:
        "docker://staphb/samtools:1.19"
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/index_sample_bam/{sample}.log",
    message:
        "Indexing bam for {wildcards.sample}"
    threads: config["threads"]["samtools"]
    resources:
        mem_gb=config["mem_gb"]["samtools"],
    shell:
        """
samtools index {input.bam} 2>&1>>{log}
        """


rule copy_sample_vcf:
    input:
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
    output:
        vcf=temp(OUT + "/mtb_typing/prepared_files/{sample}.vcf"),
    log:
        OUT + "/log/copy_sample_vcf/{sample}.log",
    shell:
        """
cp {input.vcf} {output.vcf} 2>&1>{log}
        """


rule index_sample_vcf_gatk:
    input:
        vcf=OUT + "/mtb_typing/prepared_files/{sample}.vcf",
    output:
        tbi=OUT + "/mtb_typing/prepared_files/{sample}.vcf.idx",
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/index_sample_vcf/{sample}.log",
    message:
        "Indexing vcf for {wildcards.sample}"
    threads: config["threads"]["samtools"]
    resources:
        mem_gb=config["mem_gb"]["samtools"],
    shell:
        """
gatk IndexFeatureFile -I {input.vcf} 2>&1>>{log}
        """


rule copy_ref:
    input:
        reference=lambda wildcards: SAMPLES[wildcards.sample]["reference"],
    output:
        reference=temp(OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta"),
    message:
        "Copying reference genome to output directory"
    log:
        OUT + "/log/copy_sample_bam/{sample}.log",
    shell:
        """
cp {input.reference} {output.reference}
        """


rule bwa_index_ref:
    input:
        reference=OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta",
    output:
        reference=temp(OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta.sa"),
    message:
        "Indexing reference genome for {wildcards.sample} using bwa"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/bwa:0.7.17"
    log:
        OUT + "/log/bwa_index_ref/{sample}.log",
    message:
        "Indexing ref (bwa) for {wildcards.sample}"
    threads: config["threads"]["bwa"]
    resources:
        mem_gb=config["mem_gb"]["bwa"],
    shell:
        """
bwa index {input} 2>&1>{log}
       """


rule gatk_index_ref:
    input:
        reference=OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta",
    output:
        reference=temp(OUT + "/mtb_typing/prepared_files/{sample}_ref.dict"),
    message:
        "Indexing reference genome for {wildcards.sample} using GATK"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/gatk_index_ref/{sample}.log",
    message:
        "Indexing ref (GATK) for {wildcards.sample}"
    threads: config["threads"]["gatk"]
    resources:
        mem_gb=config["mem_gb"]["gatk"],
    shell:
        """
gatk CreateSequenceDictionary -R {input.reference} 2>&1>{log}
        """


rule samtools_index_ref:
    input:
        reference=OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta",
    output:
        reference=temp(OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta.fai"),
    message:
        "Indexing reference genome for {wildcards.sample} using samtools"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.19"
    log:
        OUT + "/log/samtools_index_ref/{sample}.log",
    message:
        "Indexing ref (samtools) for {wildcards.sample}"
    threads: config["threads"]["samtools"]
    resources:
        mem_gb=config["mem_gb"]["samtools"],
    shell:
        """
samtools faidx {input.reference} 2>&1>{log}
        """


rule prepare_snpeff_config:
    input:
        template="files/mtb/snpeff_template.config",
        genbank=lambda wildcards: SAMPLES[wildcards.sample]["reference_genbank"],
    output:
        db_dir=directory(
            OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff_ref"
        ),
        config=temp(OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff.config"),
    conda:
        "../envs/biopython.yaml"
    container:
        "docker://quay.io/biocontainers/biopython:1.78"
    log:
        OUT + "/log/prepare_snpeff_config/{sample}.log",
    message:
        "Preparing SnpEff config for {wildcards.sample}"
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
cp {input.template} {output.config} 2>{log}
mkdir -p {output.db_dir} 2>>{log}
cp {input.genbank} {output.db_dir}/genes.gbk 2>>{log}
python3 workflow/scripts/prepare_snpeff.py {output.config} {input.genbank} 2>&1>>{log}
        """


rule build_snpeff_db:
    input:
        db_dir=OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff_ref",
        config=OUT + "/mtb_typing/prepared_reference_data/{sample}/snpeff.config",
    output:
        touch(OUT + "/mtb_typing/prepared_reference_data/{sample}/build_snpeff_db.done"),
    conda:
        "../envs/snpeff.yaml"
    container:
        "docker://staphb/snpeff:5.1"
    log:
        OUT + "/log/build_snpeff_db/{sample}.log",
    params:
        use_singularity=config["use_singularity"],
    message:
        "Building SnpEff db for {wildcards.sample}"
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
WORKDIR=$(dirname {input.config})
CONFIG_NAME=$(basename {input.config})
DB_NAME=$(basename {input.db_dir})
cd $WORKDIR
$EXEC build -genbank -v $DB_NAME -config $CONFIG_NAME -dataDir . 2>&1>{log} 
        """


rule prepare_ab_table:
    input:
        csv=lambda wildcards: SAMPLES[wildcards.sample]["resistance_variants_csv"],
    output:
        uncompressed=temp(
            OUT + "/mtb_typing/prepared_reference_data/{sample}/ab_table.tab"
        ),
    params:
        POS="genomepos",
        metadata=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_columns"
        ],
    log:
        OUT + "/log/prepare_ab_table/{sample}.log",
    message:
        "Coverting AMR table for {wildcards.sample}"
    shell:
        """
python workflow/scripts/convert_ab_table.py \
--force-chrom NC_000962.3 \
--POS {params.POS} \
--other {params.metadata:q} \
{input} \
{output.uncompressed} 2>&1>{log}
        """


rule compress_index_ab_table:
    input:
        uncompressed=OUT + "/mtb_typing/prepared_reference_data/{sample}/ab_table.tab",
    output:
        compressed=OUT + "/mtb_typing/prepared_reference_data/{sample}/ab_table.tab.gz",
        index=OUT + "/mtb_typing/prepared_reference_data/{sample}/ab_table.tab.gz.tbi",
    conda:
        "../envs/bcftools.yaml"
    container:
        "docker://staphb/htslib:1.17"
    log:
        OUT + "/log/compress_index_ab_table/{sample}.log",
    message:
        "Compressing and indexing AMR table for {wildcards.sample}"
    threads: config["threads"]["bcftools"]
    resources:
        mem_gb=config["mem_gb"]["bcftools"],
    shell:
        """
bgzip -c {input.uncompressed} 1> {output.compressed} 2>{log}
tabix -s 1 -b 2 -e 2 {output.compressed} 2>&1>>{log}
        """


rule generate_ab_table_header:
    output:
        OUT + "/mtb_typing/prepared_reference_data/{sample}/ab_table.header",
    params:
        columns=lambda wildcards: SAMPLES[wildcards.sample][
            "resistance_variants_columns"
        ],
    shell:
        """
python workflow/scripts/generate_ab_table_header.py \
{params.columns:q} \
{output}
        """
