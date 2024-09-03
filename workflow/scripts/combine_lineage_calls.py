# merge two dataframes and assert a single row is left

import argparse
import pandas as pd
import re

# define lineage4 sublineages from PMID 30789126
lineage4_sublineages = ["ghana", "xtype", "haarlem", "ural", "tur", "lam", "stype", "uganda", "cameroon"]

def parse_str_to_dict(lineage_str):
    lineages = lineage_str.split(",")
    parsed_data = {}
    for lineage in lineages:
        # get name of the lineage
        lineage_name = lineage.split("(")[0]
        snp_string = lineage.split("(")[1].split(")")[0]
        # get the count of identified snps and count of total snps
        snp_counts = int(snp_string.split("/")[0])
        snp_total = int(snp_string.split("/")[1])
        snp_ratio = snp_counts / snp_total
        parsed_data[lineage_name] = (snp_counts, snp_total, snp_ratio)
    return parsed_data

def get_highest_ratio(lineage_dict):
    if len(lineage_dict) == 1:
        return list(lineage_dict.keys())
    elif len(lineage_dict) == 0:
        return []
    else:
        # get the highest ratio
        max_ratio = max([lineage_dict[lineage][2] for lineage in lineage_dict.keys()])
        # get the lineage(s) with the highest ratio
        max_lineages = [lineage for lineage in lineage_dict.keys() if lineage_dict[lineage][2] == max_ratio]
        return max_lineages

def decide_type(lineage_dict, sublineages):
    output = []
    max_lineages = get_highest_ratio(lineage_dict)
    for lineage in max_lineages:
        output.append(f"{lineage}({lineage_dict[lineage][0]}/{lineage_dict[lineage][1]})")
        if lineage == "lineage4":
            # find the next highest after lineage4
            lineage_dict_copy = lineage_dict.copy()
            lineage_dict_copy.pop("lineage4")
            max_lineages_without_lineage4 = get_highest_ratio(lineage_dict_copy)
            for sublineage in max_lineages_without_lineage4:
                if sublineage in sublineages:
                    output.append(f"{sublineage}({lineage_dict[sublineage][0]}/{lineage_dict[sublineage][1]})")
    return ",".join(sorted(output))

def parse_value(value):
    """
    Parse a column from fast-lineage-caller

    Results have a format: typeA(a/b)[,typeB(x/y)]

    Where:
    - typeA and typeB are the identified lineage
    - a and x are the number of identified SNPs for lineages typeA and typeB, resp.
    - b and y are the number of lineage-specific SNPs in the schemes of typeA and typeB, resp.

    The ancestral lineage4, to which the reference genome belongs, is identified by conserved positions.
    If another type is confidently identified in addition to lineage4, this means the other type is preferred.
    Only if lineage4 is the sole confidently identified type, it should be outputted
    """
    if not isinstance(value, str):
        top_hit = value
    elif value == "NA":
        top_hit = "NA"
    else:
        parsed_lineages = parse_str_to_dict(value)
        top_hit = decide_type(parsed_lineages, lineage4_sublineages)
    return top_hit

def parse_counts(df_input):
    df = df_input.copy()
    df.columns = [f"{col}_counts" if col != "Isolate" else col for col in df.columns]
    list_warnings = []
    for col in df.columns[1:]:
        # remove _counts from column name
        col_name = col.replace("_counts", "")
        # get string value from first row, column matching col
        best_type = parse_value(df[col].iloc[0])
        df[col_name] = best_type

    return df


def main(args):
    # read in all lineage calls
    df_standard = pd.read_csv(args.standard, sep="\t")

    df_standard_parsed = parse_counts(df_standard)

    df_custom = pd.read_csv(args.custom, sep="\t")


    # merge all dataframes horizontally on Isolate
    merged_df = pd.merge(df_standard_parsed, df_custom, on="Isolate")

    assert merged_df.shape[0] == 1

    merged_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine lineage calls")

    parser.add_argument("--standard", help="Path to standard lineage calls")
    parser.add_argument("--custom", help="Path to custom lineage calls")
    parser.add_argument("--output", help="Path to output file")

    args = parser.parse_args()

    main(args)