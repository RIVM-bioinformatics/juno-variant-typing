#!/usr/bin/env python3

import io
import json

import pandas as pd


def parse_CollectWgsMetrics(path):
    lines_summary_table = []
    # Append two lines after line starting with "## METRICS CLASS" to lines_summary_table
    with open(path, "r") as f:
        for line in f:
            if line.startswith("## METRICS CLASS"):
                for _ in range(2):
                    lines_summary_table.append(next(f))
                break
    # Create a dataframe from the lines_summary_table
    df = pd.read_csv(io.StringIO("\n".join(lines_summary_table)), sep="\t")
    # Get Mean coverage
    print(df)
    mean_coverage = df["MEAN_COVERAGE"].values[0]
    median_coverage = df["MEDIAN_COVERAGE"].values[0]
    return mean_coverage, median_coverage


def parse_rrs_rrl_snp_counts(path):
    # Get value from second line of file and strip newline
    with open(path, "r") as f:
        next(f)
        rrs_rrl_snp_counts = next(f).strip()
    return rrs_rrl_snp_counts


def parse_lineage_call(path):
    # Read lineage call as df
    df = pd.read_csv(path, sep="\t")
    df = df.drop(columns=["Isolate"])
    df.fillna("-", inplace=True)
    # Get lineage calls
    lineage_dict = df.to_dict(orient="index")[0]
    return lineage_dict


def main(args):
    mean_coverage, median_coverage = parse_CollectWgsMetrics(args.picard)
    rrs_rrl_snp_counts = parse_rrs_rrl_snp_counts(args.rrs_rrl_snp_counts)
    lineage_dict = parse_lineage_call(args.lineage_call)
    # Create dictionary
    tb_dict = {
        "mean_coverage": str(mean_coverage),
        "median_coverage": str(median_coverage),
        "rrs_rrl_snp_counts": str(rrs_rrl_snp_counts),
        "lineage": lineage_dict,
    }
    # Write dictionary to JSON file
    with open(args.output, "w") as f:
        json.dump(tb_dict, f, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create a JSON file for the TB DB")

    parser.add_argument(
        "--picard", type=str, help="Picard CollectWgsMetrics output txt file"
    )
    parser.add_argument(
        "--lineage-call", type=str, help="Input rrs-rrl SNP counts file"
    )
    parser.add_argument(
        "--rrs-rrl-snp-counts", type=str, help="Input rrs-rrl SNP counts file"
    )
    parser.add_argument("--output", type=str, help="Output JSON file")

    args = parser.parse_args()

    main(args)
