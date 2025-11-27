#!/usr/bin/env python3

import argparse
import pandas as pd

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
                    row_copy["gene"] = indel_row.get("gene", "-")
                    row_copy["drug"] = indel_row.get("drug", "-")
                    row_copy["confidence"] = "INDEL gene list"
                    annotated_rows.append(row_copy)
            else:
                annotated_rows.append(row)
        else:
            annotated_rows.append(row)
    return pd.DataFrame(annotated_rows)

def main(args):
    df = pd.read_csv(args.input, sep="\t")
    df_indel_gene_list = pd.read_csv(args.reference_data, sep="\t")
    annotated_df = annotate_indels_with_gene_list(df, df_indel_gene_list)
    # Only proceed if 'confidence' column exists
    if "confidence" in annotated_df.columns:
        annotated_df = annotated_df[annotated_df["confidence"].notna() & (annotated_df["confidence"] != "")]
        if not annotated_df.empty:
            annotated_df.to_csv(args.output, sep="\t", index=False)
        else:
            # Write headers only, empty file
            annotated_df.iloc[0:0].to_csv(args.output, sep="\t", index=False)
    else:
        # Write headers only, empty file
        annotated_df.iloc[0:0].to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input raw variants tsv", required=True)
    parser.add_argument("--reference_data", help="Indel gene list tsv", required=True)
    parser.add_argument("--merge_cols", help="Comma separated columns to merge on", required=False)
    parser.add_argument("--keep_cols", help="Comma separated columns to keep from reference", required=False)
    parser.add_argument("--output", help="Output tsv", required=True)
    args = parser.parse_args()
    main(args)