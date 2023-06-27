#!/usr/bin/env python3

import pandas as pd
import argparse
import re
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("--force-chrom",
                    help="Fixed value for CHROM to be used in table",
                    default="NC_000962.3",
                    type=str)
parser.add_argument("--POS",
                    help="Column to parse genome positions",
                    type=str)
parser.add_argument("--REF",
                    help="Column to parse reference allele",
                    type=str)
parser.add_argument("--ALT",
                    help="Column to parse alternative allele",
                    type=str)
parser.add_argument("--other",
                    help="""
                    Column(s) to parse for metadata to be included.
                    If you want multiple columns, provide a comma-separated list (e.g. --other ab,ab_class)
                    """,
                    type=str)
parser.add_argument(help="Input csv file",
                    dest="input",
                    type=Path)
parser.add_argument(help="Output tab file",
                    dest="output",
                    type=Path)

args = parser.parse_args()

df_in = pd.read_csv(args.input, sep=';')

metadata_cols = args.other.split(',')

ordered_columns = [args.POS] + metadata_cols

for col in ordered_columns:
    assert col in df_in.columns, f"Column {col} cannot be found in {args.input}"

df_out = df_in[ordered_columns]

df_out.insert(0, '#CHROM', args.force_chrom)

df_out = df_out.sort_values('genomepos')
df_out = df_out.fillna("missing")

df_out.rename(columns = {'genomepos': 'POS'}, inplace = True)

df_out.to_csv(args.output, sep='\t', index=False)