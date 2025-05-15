#!/usr/bin/env python3
import pandas as pd
import re
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Process Iso-Seq transcript model files.")
parser.add_argument("--gtf", required=True, help="Path to the input GTF file")
parser.add_argument("--tsv", required=True, help="Path to the grouped counts TSV file")
parser.add_argument("--output", default="06_isoquant.flnc.count.csv", help="Output CSV file name")
args = parser.parse_args()

gtf_in = Path(args.gtf)
tsv_file = Path(args.tsv)
csv_file = Path(args.output)

# Step 1: Replace header in GTF
with gtf_in.open("r") as f:
    gtf_lines = [re.sub(r"^# ", "##", line) for line in f]
with gtf_in.open("w") as f:
    f.writelines(gtf_lines)

# Step 2: Convert TSV to CSV
with tsv_file.open("r") as tsv, csv_file.open("w") as csv:
    for line in tsv:
        line = line.replace("#feature_id", "id")
        csv.write(line.replace("\t", ","))

# Step 3: Extract transcript IDs from GTF
with gtf_in.open("r") as f:
    gtf_ids = set(
        re.search(r'transcript_id "([^"]+)"', line).group(1)
        for line in f
        if not line.startswith("#") and "\ttranscript\t" in line
    )

# Step 4: Extract transcript IDs from CSV (in memory)
df = pd.read_csv(csv_file)
flnc_ids = set(df['id'].dropna().astype(str))

# Step 5: Find missing transcript IDs
missing_ids = sorted(gtf_ids - flnc_ids)

# Step 6: Add missing transcripts with 0 counts
num_cols = len(df.columns)
with csv_file.open("a") as out:
    for tid in missing_ids:
        out.write(tid + "," + ",".join(["0"] * (num_cols - 1)) + "\n")
