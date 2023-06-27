#!/usr/bin/env python3

import pathlib
from Bio import SeqIO

def check_codontable(codontable, config):
    with open(config, "r") as file:
        available_codontables = []
        for line in file:
            if line.startswith("codon."):
                codontable_name = line.split('\t')[0].split('.')[1]
                available_codontables.append(codontable_name)
    assert codontable in available_codontables, f"Codon table {codontable} is not available in snpEff template config {config}"

def get_chromosomes(genbank_ref_path):
    list_chromosomes = []
    with open(genbank_ref_path, "r") as genbank_ref:
        for record in SeqIO.parse(genbank_ref, "genbank"):
            list_chromosomes.append(record.id)
    return list_chromosomes


def main(args):
    check_codontable(args.codontable, args.config)
    list_chromosomes = get_chromosomes(args.genbank)
    comma_str_chromosomes = ", ".join(list_chromosomes)
    with open(args.config, "a+") as file:
        file.write("\n")
        file.write("# Reference\n")
        file.write("snpeff_ref.genome : SnpEff reference\n")
        file.write(f"snpeff_ref.chromosomes : {comma_str_chromosomes}\n")
        for chrom in list_chromosomes:
            file.write(f"snpeff_ref.{chrom}.codonTable : {args.codontable}\n")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument(help="SnpEff config file",
                        dest="config",
                        metavar="SNPEFF_CONFIG",
                        type=pathlib.Path)
    parser.add_argument(help="Genbank reference file",
                        dest="genbank",
                        metavar="GENBANK_REF",
                        type=pathlib.Path)
    parser.add_argument("--codontable",
                        help="Codon table to use. Should be present in SnpEff template config [Bacterial_and_Plant_Plastid]",
                        type=str,
                        metavar="STR",
                        default="Bacterial_and_Plant_Plastid")

    args = parser.parse_args()

    main(args)