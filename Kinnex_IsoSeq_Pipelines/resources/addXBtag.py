#!/usr/bin/env python3

import pysam
import argparse
from pathlib import Path

def process_bam(input_bam, sample_name, output_bam):
    infile = pysam.AlignmentFile(input_bam, "rb", check_sq=False)

    # Modify @RG headers to set SM:sample_name
    header = infile.header.to_dict()
    if "RG" in header:
        for rg in header["RG"]:
            rg["SM"] = sample_name

    outfile = pysam.AlignmentFile(output_bam, "wb", header=header)

    prefix = f"{sample_name}_"
    for read in infile.fetch(until_eof=True):
        read.query_name = prefix + read.query_name

        try:
            cb = read.get_tag("CB")
            read.set_tag("XB", prefix + cb, value_type='Z')
        except KeyError:
            pass

        outfile.write(read)

    infile.close()
    outfile.close()

def main():
    parser = argparse.ArgumentParser(description="Prefix read names, copy CB to XB, and update SM tag in header")
    parser.add_argument("input_bam", help="Input BAM")
    parser.add_argument("sample_name", help="Sample name")
    parser.add_argument("output_bam", help="Output BAM (default: <input>.tagged.bam)", default=None)
    args = parser.parse_args()

    input_path = Path(args.input_bam)
    output_path = args.output_bam or input_path.with_suffix(".tagged.bam")

    process_bam(str(input_path), args.sample_name, str(output_path))
    print(f"Output written to: {output_path}")

if __name__ == "__main__":
    main()
