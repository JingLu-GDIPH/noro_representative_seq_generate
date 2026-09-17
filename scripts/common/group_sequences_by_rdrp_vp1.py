#!/usr/bin/env python3
"""Group norovirus sequences by identical RdRp and VP1 genotype labels.

Input FASTA identifiers are expected to start with:

    <RdRp genotype>_<VP1 genotype>_<accession>_<collection date>

For example:

    GII.P16_GII.4_AB123456_2012

The script produces files compatible with the existing downstream pipeline:
``cluster_info.csv`` and ``clustered_sequences.fasta``.
"""

import argparse
import csv
import re
from collections import defaultdict

from Bio import SeqIO


RDRP_PATTERN = re.compile(r"^G(?:I|II|IX)\.P[A-Za-z0-9]+$")
VP1_PATTERN = re.compile(r"^G(?:I|II|IX)\.[A-Za-z0-9]+$")


def parse_genotype_pair(sequence_id):
    """Return (rdrp, vp1) parsed from the first two FASTA ID fields."""
    parts = sequence_id.split("_")
    if len(parts) < 2:
        return "UnknownRdRp", "UnknownVP1"

    rdrp = parts[0].strip() or "UnknownRdRp"
    vp1 = parts[1].strip() or "UnknownVP1"

    if not RDRP_PATTERN.match(rdrp):
        rdrp = "UnknownRdRp"
    if not VP1_PATTERN.match(vp1):
        vp1 = "UnknownVP1"

    return rdrp, vp1


def group_records(input_fasta, cluster_info, sequences_output, report_output):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    groups = defaultdict(list)

    for record in records:
        rdrp, vp1 = parse_genotype_pair(record.id)
        groups[(rdrp, vp1)].append(record)

    ordered_groups = sorted(
        groups.items(),
        key=lambda item: (-len(item[1]), item[0][0], item[0][1]),
    )

    rows = []
    report_rows = []
    for index, ((rdrp, vp1), group_records_) in enumerate(ordered_groups, start=1):
        cluster_id = f"cluster_{index:04d}"
        genotype_pair = f"{rdrp}_{vp1}"
        report_rows.append(
            {
                "cluster_id": cluster_id,
                "rdrp_genotype": rdrp,
                "vp1_genotype": vp1,
                "genotype_pair": genotype_pair,
                "sequence_count": len(group_records_),
            }
        )

        for record in group_records_:
            rows.append(
                {
                    "sequence_id": record.id,
                    "cluster_id": cluster_id,
                    "rdrp_genotype": rdrp,
                    "vp1_genotype": vp1,
                    "genotype_pair": genotype_pair,
                    "sequence_length": len(str(record.seq).replace("-", "")),
                }
            )

    SeqIO.write(records, sequences_output, "fasta")

    with open(cluster_info, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sequence_id",
                "cluster_id",
                "rdrp_genotype",
                "vp1_genotype",
                "genotype_pair",
                "sequence_length",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    with open(report_output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "cluster_id",
                "rdrp_genotype",
                "vp1_genotype",
                "genotype_pair",
                "sequence_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(report_rows)

    print(f"Grouped {len(records)} records into {len(ordered_groups)} RdRp+VP1 clusters")
    print(f"Cluster info: {cluster_info}")
    print(f"Cluster report: {report_output}")


def main():
    parser = argparse.ArgumentParser(
        description="Group FASTA records by identical RdRp and VP1 genotype labels."
    )
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--cluster_info", required=True, help="Output cluster_info.csv")
    parser.add_argument("--sequences_output", required=True, help="Output FASTA")
    parser.add_argument("--report", required=True, help="Output cluster summary TSV")
    args = parser.parse_args()

    group_records(args.input, args.cluster_info, args.sequences_output, args.report)


if __name__ == "__main__":
    main()
