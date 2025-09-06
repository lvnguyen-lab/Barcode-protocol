#!/usr/bin/env python

import sys
import os
import pandas as pd
from itertools import combinations

indir = sys.argv[1]
outdir = sys.argv[2]

files = os.listdir(indir)

# Set up list to hold dictionaries with R1 and R2 files
paired_files = []

# Pair files
all_combos = list(combinations(files, 2))

for pair in all_combos:
    pair_1_split = pair[0].split("_")
    pair_2_split = pair[1].split("_")
    sample = "_".join(pair_1_split[0:3])
    if pair_1_split[0:3] == pair_2_split[0:3]:
        if pair_1_split[3] == "R1":
            paired_dict = {
                "sample": sample,
                "R1": pair[0],
                "R2": pair[1]
            }
        else:
            paired_dict = {
                "sample": sample,
                "R1": pair[1],
                "R2": pair[0]                
            }
        paired_files.append(paired_dict)

paired_df = pd.DataFrame(paired_files)

paired_df.to_csv(f"{outdir}/paired_fastq.csv", index=False, header=False)