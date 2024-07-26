# merge two dataframes and assert a single row is left

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Combine lineage calls")

parser.add_argument("lineage_calls", help="Path to lineage calls", nargs="+")
parser.add_argument("--output", help="Path to output file")

args = parser.parse_args()

# read in all lineage calls
dfs = []
for lineage_call in args.lineage_calls:
    df = pd.read_csv(lineage_call, sep="\t")
    dfs.append(df)

# merge all dataframes horizontally on Isolate
merged_df = dfs[0]
for df in dfs[1:]:
    merged_df = pd.merge(merged_df, df, on="Isolate")

assert merged_df.shape[0] == 1

merged_df.to_csv(args.output, sep="\t", index=False)
