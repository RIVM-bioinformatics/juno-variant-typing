def choose_species(wildcards):
    if (SAMPLES[wildcards.sample]["genus"] == "mycobacterium") & (
        SAMPLES[wildcards.sample]["species"] == "tuberculosis"
    ):
        return [
            OUT + "/mtb_typing/lineage_call/{sample}.tsv",
            OUT + "/mtb_typing/contamination_check/coll_positions/{sample}.tsv",
            OUT + "/mtb_typing/contamination_check/rrs_rrl_contamination/{sample}.tsv",
        ]
    else:
        return OUT + "/typing_check/{sample}/no_typing_necessary.txt"


rule aggregate_species:
    input:
        choose_species,
    output:
        temp(OUT + "/typing_check/{sample}_done.txt"),
    message:
        "Checking correct typing ran for {wildcards.sample}"
    threads: 1
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        "touch {output}"


rule no_typing:
    input:
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
    output:
        temp(OUT + "/typing_check/{sample}/no_typing_necessary.txt"),
    message:
        "Skipping typing step for {wildcards.sample}."
    threads: 1
    shell:
        """
        touch {output}
        """
