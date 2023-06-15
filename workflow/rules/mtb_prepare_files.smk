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
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/index_sample_bam/{sample}.log",
    shell:
        """
samtools index {input.bam} 2>&1>>{log}
        """


rule copy_ref:
    input:
        reference=lambda wildcards: SAMPLES[wildcards.sample]["reference"],
    output:
        reference=temp(OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta"),
    message:
        "Copying reference genome to output directory"
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
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
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
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
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
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/samtools_index_ref/{sample}.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
samtools faidx {input.reference} 2>&1>{log}
        """
