#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Input tsv file", type=Path)
parser.add_argument("--reference_data", help="Reference data csv file", type=Path)
parser.add_argument(
    "--merge_cols", help="Comma separated list of columns to merge on", type=str
)
parser.add_argument(
    "--keep_cols",
    help="Comma separated list of columns to keep from reference data",
    type=str,
)
parser.add_argument("--output", help="Output tsv file", type=Path)
parser.add_argument("--indel-gene-list", help="Indel gene list tsv file", type=Path)
args = parser.parse_args()

rename_dict = {
    "variant_common_name": "reslist_variant",
}


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

def annotate_indels_with_gene_list(df, df_indel_gene_list):
    """
    Annotate INDEL variants with gene and drug info from indel gene list.
    """
    annotated_rows = []
    for _, row in df.iterrows():
        if row.get("TYPE") == "INDEL":
            matches = df_indel_gene_list[
                (row["POS"] >= df_indel_gene_list["start"]) &
                (row["POS"] <= df_indel_gene_list["end"])
            ]
            if not matches.empty:
                for _, indel_row in matches.iterrows():
                    row_copy = row.copy()
                    row_copy["gene"] = indel_row["gene"]
                    row_copy["drug"] = indel_row["drug"]
                    row_copy["confidence"] = "INDEL gene list"
                    annotated_rows.append(row_copy)
            else:
                annotated_rows.append(row)
        else:
            annotated_rows.append(row)
    return pd.DataFrame(annotated_rows)


def main(args):
    df = pd.read_csv(args.input, sep="\t")
    df_ref = pd.read_csv(args.reference_data, sep=",")
    # Keep only the columns that are needed
    all_keep_cols = args.merge_cols.split(",") + args.keep_cols.split(",")
    df_merged = pd.merge(
        df, df_ref[all_keep_cols], on=args.merge_cols.split(","), how="left"
    )
    df_merged_eff_parsed = parse_eff_field(df_merged)
    df_merged_eff_parsed_renamed = rename_columns(df_merged_eff_parsed, rename_dict)
    df_final = df_merged_eff_parsed_renamed.fillna("-").replace("", "-")

    #Annotate INDELs if indel gene list is provided
    if args.indel_gene_list is not None:
        df_indel_gene_list = pd.read_csv(args.indel_gene_list, sep="\t")
        df_final = annotate_indels_with_gene_list(df_final, df_indel_gene_list)
    df_final.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main(args)
