# merge two dataframes and assert a single row is left

import argparse
import pandas as pd
import re


def parse_value(value):
    if not isinstance(value, str):
        outcome = value
    elif value == "NA":
        outcome = "NA"
    elif len(value.split(",")) == 1:
        # remove string matching pattern "([0-9]+/[0-9]+)" from value
        outcome = re.sub(r"\([0-9]+/[0-9]+\)", "", value)
    else:
        # make a list of identified types
        list_of_types = value.split(",")
        
        max_ratio = -1
        list_types_with_max_ratio = []

        for item in list_of_types:
            # get type name by removing count matching pattern "([0-9]+/[0-9]+)"
            type_name = re.sub(r"\([0-9]+/[0-9]+\)", "", item)
            # get count matching pattern "([0-9]+/[0-9]+)"
            count_str = re.search(r"\([0-9]+/[0-9]+\)", item).group()
            # get first and second number from count_str,
            # indicating identified nr of SNPs and total nr of SNPs for that type
            counts = count_str.strip("()").split("/")
            ratio_type = int(counts[0]) / int(counts[1])
            if ratio_type > max_ratio:
                max_ratio = ratio_type
                list_types_with_max_ratio = [type_name]
            elif ratio_type == max_ratio:
                list_types_with_max_ratio.append(type_name)

        # get type with highest ratio
        outcome = ",".join(list_types_with_max_ratio)

    return outcome

def parse_counts(df_input):
    df = df_input.copy()
    df.columns = [f"{col}_counts" if col != "Isolate" else col for col in df.columns]
    for col in df.columns[1:]:
        # remove _counts from column name
        col_name = col.replace("_counts", "")
        # get string value from first row, column matching col
        df[col_name] = parse_value(df[col].iloc[0])
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