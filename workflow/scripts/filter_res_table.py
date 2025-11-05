#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=Path)
    parser.add_argument("--ab-column", type=str)
    parser.add_argument("--indel-gene-list", type=Path)
    parser.add_argument("--output", type=Path)

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str, na_values=["-"])
    df["POS"] = df["POS"].astype(int)
    df_indel_gene_list = pd.read_csv(args.indel_gene_list, sep="\t")

    # just loop through all variants for script clarity
    list_rows_to_keep = []
    for _variant_index, variant_row in df.iterrows():
        # ab_column is filled from the WHO catalogue so these are always kept
        if pd.notnull(variant_row[args.ab_column]):
            list_rows_to_keep.append(variant_row)
        # if not detected using WHO catalogue, check whether the variant type is INDEL
        elif variant_row["TYPE"] == "INDEL":
            # and then check whether the indel variant is within start & end of listed genes
            df_indel_genes_matching_variant_row = df_indel_gene_list[
                (variant_row["POS"] >= df_indel_gene_list["start"])
                & (variant_row["POS"] <= df_indel_gene_list["end"])
            ]
            # for every matching gene in the indel gene list, output the matching variant with gene and drug values taken from the indel gene list
            # this can be useful if a variant matches two overlapping genes which may cause resistance to different drugs (theoretical edge case)
            for (
                _indel_index,
                indel_row,
            ) in df_indel_genes_matching_variant_row.iterrows():
                variant_row_copy = variant_row.copy()
                variant_row_copy["gene"] = indel_row["gene"]
                variant_row_copy["drug"] = indel_row["drug"]
                variant_row_copy["confidence"] = "INDEL gene list"
                list_rows_to_keep.append(variant_row_copy)

    df_out = pd.DataFrame(list_rows_to_keep)
    df_out.to_csv(args.output, sep="\t", index=False)
