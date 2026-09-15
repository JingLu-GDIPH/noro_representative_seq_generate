#!/usr/bin/env python3
"""Keep only vclust pairwise comparisons within the same VP1 genotype."""

import argparse
import csv
import re


GENOTYPE_PATTERN = re.compile(r"(?<![A-Z0-9])G(?:I|II|IX)\.\d+(?!\d)")


def extract_norovirus_genotype(text):
    """Extract GI.x, GII.x, or GIX.x from a sequence identifier."""
    match = GENOTYPE_PATTERN.search(str(text))
    return match.group() if match else None


def filter_pair_table(input_path, output_path):
    """Write only typed pairs whose query and reference genotypes match."""
    kept = 0
    removed = 0

    with open(input_path, newline="") as input_handle, open(
        output_path, "w", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if not reader.fieldnames or "query" not in reader.fieldnames or "reference" not in reader.fieldnames:
            raise ValueError("vclust table must contain query and reference columns")

        writer = csv.DictWriter(
            output_handle,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for row in reader:
            query_genotype = extract_norovirus_genotype(row["query"])
            reference_genotype = extract_norovirus_genotype(row["reference"])
            if query_genotype and query_genotype == reference_genotype:
                writer.writerow(row)
                kept += 1
            else:
                removed += 1

    return {"kept": kept, "removed": removed}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    stats = filter_pair_table(args.input, args.output)
    print(f"Kept same-genotype pairs: {stats['kept']}")
    print(f"Removed cross-genotype or untyped pairs: {stats['removed']}")


if __name__ == "__main__":
    main()
