#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

def read_lines(file_path):
    """Read non-empty, stripped lines from a file."""
    return [line.strip() for line in file_path.read_text().splitlines() if line.strip()]

def quote_list(items):
    """Convert a Python list of strings to a YAML-style quoted list."""
    return "[\n      " + ",\n      ".join(f'"{item}"' for item in items) + "\n    ]"

def main():
    parser = argparse.ArgumentParser(description="Generate YAML-like file from BAM paths and labels.")
    parser.add_argument("-b", "--bam_list", required=True, type=Path, help="Text file with BAM file paths.")
    parser.add_argument("-l", "--labels", required=True, type=Path, help="Text file with corresponding labels.")
    parser.add_argument("-o", "--output", default="output.yaml", type=Path, help="Output YAML file name.")
    parser.add_argument("-e", "--experiment-name", default="Experiment1", help="Name of the experiment.")

    args = parser.parse_args()

    bam_files = read_lines(args.bam_list)
    labels = read_lines(args.labels)

    if len(bam_files) != len(labels):
        sys.exit(f"ERROR: Number of BAM files ({len(bam_files)}) does not match number of labels ({len(labels)}).")

    content = f"""[
  data format: "bam",
  {{
    name: "{args.experiment_name}",
    long read files: {quote_list(bam_files)},
    labels: {quote_list(labels)},
  }}
]
"""

    args.output.write_text(content)
    print(f"✅ YAML-like file written to: {args.output}")

if __name__ == "__main__":
    main()
