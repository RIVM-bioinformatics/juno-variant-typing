#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def process_df(input_path: str, ab_column: str, output_path: str) -> None:
    """
    Filter out table rows if these have NA/NaN in both RESISTANCE_AB and RESISTANCE_AB_CLASS columns
    """
    df = pd.read_csv(input_path, sep="\t", dtype=str, na_values=["-"])
    filter_ab = ~df[ab_column].isna()
    df = df[filter_ab]
    df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=Path)
    parser.add_argument("--ab-column", type=str)
    parser.add_argument("--output", type=Path)

    args = parser.parse_args()

    process_df(args.input, args.ab_column, args.output)
