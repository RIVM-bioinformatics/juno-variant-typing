rule audit_version_gatk:
    output:
        OUT + "/version_audit/gatk.txt",
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    conda:
        "../envs/gatk_picard.yaml"
    shell:
        "gatk --version > {output}"


rule audit_version_biopython:
    output:
        OUT + "/version_audit/biopython.txt",
    container:
        "docker://quay.io/biocontainers/biopython:1.78"
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python3 -c 'import Bio; print(Bio.__version__)' > {output}
        """


rule audit_version_bcftools:
    output:
        OUT + "/version_audit/bcftools.txt",
    container:
        "docker://staphb/bcftools:1.19"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools --version > {output}"


rule audit_version_bedtools:
    output:
        OUT + "/version_audit/bedtools.txt",
    container:
        "docker://staphb/bedtools:2.30.0"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bedtools --version > {output}"


rule audit_version_bwa:
    output:
        OUT + "/version_audit/bwa.txt",
    container:
        "docker://staphb/bwa:0.7.17"
    conda:
        "../envs/bwa_samtools.yaml"
    shell:
        "bwa 2>&1 | grep Version > {output} || true"


rule audit_version_bgzip:
    output:
        OUT + "/version_audit/bgzip.txt",
    container:
        "docker://staphb/htslib:1.17"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bgzip --version > {output}"


rule audit_version_samtools:
    output:
        OUT + "/version_audit/samtools.txt",
    container:
        "docker://staphb/samtools:1.19"
    conda:
        "../envs/bwa_samtools.yaml"
    shell:
        "samtools --version > {output}"


rule audit_version_seqkit:
    output:
        OUT + "/version_audit/seqkit.txt",
    container:
        "docker://staphb/seqkit:2.7.0"
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit version > {output}"


rule audit_version_snpeff:
    output:
        OUT + "/version_audit/snpeff.txt",
    container:
        "docker://staphb/snpeff:5.1"
    conda:
        "../envs/snpeff.yaml"
    params:
        use_singularity=config["use_singularity"],
    shell:
        """
if [ {params.use_singularity} == "True" ]
then
    EXEC=snpeff
else
    EXEC=snpEff
fi
$EXEC -version > {output}
        """


rule combine_version:
    input:
        OUT + "/version_audit/gatk.txt",
        OUT + "/version_audit/biopython.txt",
        OUT + "/version_audit/bcftools.txt",
        OUT + "/version_audit/bedtools.txt",
        OUT + "/version_audit/bwa.txt",
        OUT + "/version_audit/bgzip.txt",
        OUT + "/version_audit/samtools.txt",
        OUT + "/version_audit/seqkit.txt",
        OUT + "/version_audit/snpeff.txt",
    output:
        OUT + "/version_audit/versions.txt",
    shell:
        "cat {input} > {output}"
