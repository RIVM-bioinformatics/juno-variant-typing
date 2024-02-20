import yaml


sample_sheet = config["sample_sheet"]
with open(sample_sheet) as f:
    SAMPLES = yaml.safe_load(f)

for param in ["threads", "mem_gb"]:
    for k in config[param]:
        config[param][k] = int(config[param][k])

# print(SAMPLES)

OUT = config["output_dir"]
INPUT = config["input_dir"]


localrules:
    all,
    no_typing,
    aggregate_species,
    copy_sample_bam,
    copy_ref,
    prepare_ab_table,
    generate_ab_table_header,
    mtb_filter_res_table_positions,


include: "workflow/rules/choose_species.smk"
include: "workflow/rules/mtb_prepare_files.smk"
include: "workflow/rules/mtb_typing.smk"


rule all:
    input:
        expand(OUT + "/typing_check/{sample}_done.txt", sample=SAMPLES),
        expand(
            OUT + "/mtb_typing/annotated_resistance_filtered/{sample}.tsv",
            sample=SAMPLES,
        ),
