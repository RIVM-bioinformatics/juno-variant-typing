#!/usr/bin/env python3

import sys

metadata = sys.argv[1]

command = []

for item in metadata.split(","):
    command.append(f"-F {item}")

command_cat = " ".join(command)

print(command_cat)
