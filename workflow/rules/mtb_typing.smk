rule mtb_lineage_id:
    input:
        vcf = lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
    output:
        tsv = OUT + "/mtb_typing/lineage_call/{sample}.tsv",
    conda:
        "../envs/fast_lineage_caller.yaml"
    log:
        OUT + "/log/lineage_id/{sample}.log"
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
        bam = OUT + "/mtb_typing/prepared_files/{sample}.bam",
        bai = OUT + "/mtb_typing/prepared_files/{sample}.bam.bai",
        reference = OUT + "/mtb_typing/prepared_files/{sample}_ref.fasta",
        dummy = OUT + "/mtb_typing/prepared_files/{sample}_ref.dict",
        bed = lambda wildcards: SAMPLES[wildcards.sample]["single_copy_bed"],
    output:
        tsv = OUT + "/mtb_typing/contamination_check/coll_positions/{sample}.tsv"
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/mtb_typing/coll_contamination/{sample}.log"
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
        vcf = lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
    output:
        tsv = OUT + "/mtb_typing/contamination_check/rrs_rrl_contamination/{sample}.tsv",
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/mtb_typing/rrs_rrl_contamination/{sample}.log"
    shell:
        """
gatk CountVariants \
-V {input.vcf} \
1> {output.tsv} \
2> {log}
        """