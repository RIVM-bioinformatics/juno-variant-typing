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

