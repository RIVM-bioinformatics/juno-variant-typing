import yaml

sample_sheet = config["sample_sheet"]
with open(sample_sheet) as f:
    SAMPLES = yaml.safe_load(f)

OUT = config["output_dir"]


wildcard_constraints:
    sample="[^/]{1,100}",


rule mark_variants_by_proximity:
    input:
        vcf="",
    output:
        bed="",
    conda:
        "../envs/bcftools.yaml"
    container:
        "docker://staphb/bedtools:2.30.0"
    log:
        "",
    params:
        proximity_threshold=lambda wildcards: SAMPLES[wildcards.sample]["consensus"][
            "proximity_threshold"
        ],
    message:
        "Marking regions with variants too close together for {wildcards.sample}"
    threads: config["threads"]["bedtools"]
    resources:
        mem_gb=config["mem_gb"]["bedtools"],
    shell:
        """
# merged VCF features with more than one entry within x bp are written to a bed file for filtering
bedtools merge -d {params.proximity_threshold} -c 1 -o count -i {input.vcf} 2>{log} |\
awk '$4 > 1 {{print $0}}' 1>{output.bed} 2>>{log}
        """


rule subset_fixed_snps_mnps_from_vcf:
    input:
        vcf="",
    output:
        vcf="",
    conda:
        "../envs/bcftools.yaml"
    container:
        "docker://staphb/bcftools:1.19"
    log:
        "",
    params:
        AF_threshold=lambda wildcards: SAMPLES[wildcards.sample]["consensus"][
            "AF_threshold"
        ],
    message:
        "Subsetting fixed snps from vcf for {wildcards.sample}"
    threads: config["threads"]["bcftools"]
    resources:
        mem_gb=config["mem_gb"]["bcftools"],
    shell:
        """
bcftools filter --include "FORMAT/AF>={params.AF_threshold}" {input.vcf} 2>{log} |\
bcftools filter --include "TYPE='snp' | TYPE='mnp'" \
1>{output.vcf} \
2>>{log}
        """


rule subset_low_confidence_variants_from_vcf:
    input:
        vcf="",
    output:
        vcf="",
    conda:
        "../envs/bcftools.yaml"
    container:
        "docker://staphb/bcftools:1.19"
    log:
        "",
    params:
        AF_threshold=lambda wildcards: SAMPLES[wildcards.sample]["consensus"][
            "AF_threshold"
        ],
        tlod_threshold=lambda wildcards: SAMPLES[wildcards.sample]["consensus"][
            "tlod_threshold"
        ],
    message:
        "Subsetting unfixed snps from vcf for {wildcards.sample}"
    threads: config["threads"]["bcftools"]
    resources:
        mem_gb=config["mem_gb"]["bcftools"],
    shell:
        """
bcftools filter \
--include "FORMAT/AF < {params.AF_threshold} | TLOD < {params.tlod_threshold}" \
{input.vcf} \
1>{output.vcf} \
2>{log}
        """


rule zip_and_index_sample_vcf_bcftools:
    input:
        vcf="",
    output:
        vcf_gz="",
        tbi="",
    conda:
        "../envs/bcftools.yaml"
    container:
        "docker://staphb/htslib:1.17"
    log:
        "",
    message:
        "Compressing and indexing vcf for {wildcards.sample}"
    threads: config["threads"]["bcftools"]
    resources:
        mem_gb=config["mem_gb"]["bcftools"],
    shell:
        """
bgzip -c {input.vcf} 1> {output.vcf_gz} 2>{log}
tabix -p vcf {output.vcf_gz} 2>&1>>{log}
        """


rule introduce_mutations_to_reference:
    input:
        vcf_gz="",
        tbi="",
        reference="",
    output:
        fasta="",
    conda:
        "../envs/bcftools.yaml"
    container:
        "docker://staphb/bcftools:1.19"
    log:
        "",
    message:
        "Creating consensus for {wildcards.sample}"
    threads: config["threads"]["bcftools"]
    resources:
        mem_gb=config["mem_gb"]["bcftools"],
    shell:
        """
bcftools consensus \
-f {input.reference} \
-o {output.fasta} \
-s - \
{input.vcf_gz} \
2>&1>{log}
        """


rule mask_fasta_on_depth_from_bam:
    input:
        bam="",
        bai="",
        fasta="",
    output:
        fasta="",
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://staphb/bedtools:2.30.0"
    log:
        "",
    message:
        "Masking low depth regions of {wildcards.sample} consensus genome"
    params:
        masking_depth=lambda wildcards: SAMPLES[wildcards.sample]["consensus"][
            "masking_depth"
        ],
    threads: config["threads"]["bedtools"]
    resources:
        mem_gb=config["mem_gb"]["bedtools"],
    shell:
        """
bedtools genomecov -ibam {input.bam} -bga 2> {log} |\
awk -v DEPTH={params.masking_depth} '$4 < DEPTH {{print $0}}' 2>>{log} |\
bedtools maskfasta -fi {input.fasta} -bed stdin -fo {output.fasta} 2>> {log}
        """


rule mask_fasta_based_on_bed_or_vcf:
    input:
        features="",
        fasta="",
    output:
        fasta="",
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://staphb/bedtools:2.30.0"
    log:
        "",
    message:
        "Masking heterozygous variants in {wildcards.sample}"
    threads: config["threads"]["bedtools"]
    resources:
        mem_gb=config["mem_gb"]["bedtools"],
    shell:
        """
bedtools maskfasta \
-fi {input.fasta} \
-bed {input.features} \
-fo {output.fasta} 2>&1>{log}
        """


rule replace_fasta_header:
    input:
        fasta="",
    output:
        fasta="",
    conda:
        "../envs/seqkit.yaml"
    container:
        "docker://staphb/seqkit:2.7.0"
    log:
        "",
    message:
        "Replacing fasta header for {wildcards.sample}"
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
seqkit replace \
-p ^ \
-r '{wildcards.sample}_contig{{nr}} ' \
{input.fasta} \
1> {output.fasta} \
2>{log}
        """
