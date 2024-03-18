#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Input tsv file", type=Path)
parser.add_argument("--output", help="Output tsv file", type=Path)
args = parser.parse_args()

rename_dict = {
    "variant_common_name": "reslist_variant",
}


def read_df(path):
    df = pd.read_csv(path, sep="\t")
    return df


def parse_eff_field(df):
    """
    Parse the EFF field into a more readable format

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe

    Returns
    -------
    pandas.DataFrame
        Dataframe with parsed EFF field
    """
    eff_cols = [
        "impact",
        "functional_class",
        "detected_codon_change",
        "detected_amino_acid_change",
        "ref_genome_gene_name",
        "biotype",
        "gene_coding",
        "locus_tag",
        "exon_rank",
        "genotype_number",
        "warnings",
        "errors",
    ]  # variable number of columns: https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files

    unneeded_columns = [
        "EFF",
        "gene_coding",
        "exon_rank",
        "genotype_number",
        "warnings",
        "errors",
    ]

    # Get the mutation type from the EFF field
    # The mutation type is in all caps with possibly underscores before '('
    df["mutation_type"] = df["EFF"].str.extract(r"([A-Z_]+)\(")
    # extract snpeff fields between parentheses
    df_snpeff_fields = df["EFF"].str.extract("\(([^)]+)\)")
    # split fields on pipe character
    df_snpeff_fields_split = df_snpeff_fields[0].str.split("|", expand=True)
    # rename columns
    n_cols = len(df_snpeff_fields_split.columns)
    df_snpeff_fields_split.columns = eff_cols[:n_cols]
    # join back to original dataframe
    df = df.join(df_snpeff_fields_split)
    # remove unneeded columns
    for col in unneeded_columns:
        if col in df.columns:
            if (col == "warnings") or (col == "errors"):
                logging.warning("Column {col} was dropped from snpeff output")
            df = df.drop(col, axis=1)

    return df


def rename_columns(df, rename_dict):
    """
    Rename columns to more descriptive names

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe

    Returns
    -------
    pandas.DataFrame
        Renamed dataframe
    """
    df = df.rename(columns=rename_dict)
    # Check if any column ends with .AF, and if so rename the whole column name to AF
    # throw error if there are multiple columns ending with .AF
    af_cols = [col for col in df.columns if col.endswith(".AF")]
    if len(af_cols) > 1:
        raise ValueError(f"Multiple columns ending with .AF: {af_cols}")
    elif len(af_cols) == 1:
        df = df.rename(columns={af_cols[0]: "AF"})

    return df


def main(args):
    df = read_df(args.input)
    df_eff_parsed = parse_eff_field(df)
    df_eff_parsed_renamed = rename_columns(df_eff_parsed, rename_dict)
    df_final = df_eff_parsed_renamed.fillna("-").replace("", "-")
    df_final.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main(args)
