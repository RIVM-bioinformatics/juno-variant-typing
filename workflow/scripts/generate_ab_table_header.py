#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("columns",
                    type=str)
parser.add_argument("output",
                    type=str)

args = parser.parse_args()

list_cols = args.columns.split(',')

with open(args.output, "w+") as file:
    for col in list_cols:
        file.write(f"##INFO=<ID={col},Number=1,Type=String,Description=\"{col}\">\n")