#!/usr/bin/env python3
"""Rename downloaded FASTA records from VP1-typed metadata."""

import argparse
import csv
import os
import re
import tempfile
from pathlib import Path


def safe_collection_date(value):
    value = re.sub(r"[^A-Za-z0-9]+", "-", value.strip())
    return value.strip("-")


def build_seqname(row):
    genotype = row.get("vp1_genotype", "").strip()
    if not genotype:
        genotype = f"{row['group']}.Unassigned"
    date = safe_collection_date(row["collection_date"])
    return f"{genotype}_{row['accession']}_{date}"


def accession_from_fasta_id(identifier):
    match = re.search(
        r"(?:^|_)([A-Z]{1,4}_?\d{5,9})(?:\.\d+)?(?=_|$)",
        identifier,
    )
    if not match:
        raise ValueError(f"Cannot extract accession from FASTA ID: {identifier}")
    return match.group(1)


def read_metadata(path):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        columns = list(reader.fieldnames)

    if "seqname" not in columns:
        accession_index = columns.index("accession")
        columns.insert(accession_index + 1, "seqname")

    seqnames = set()
    by_accession = {}
    for row in rows:
        seqname = build_seqname(row)
        if seqname in seqnames:
            raise ValueError(f"Duplicate seqname generated: {seqname}")
        seqnames.add(seqname)
        row["seqname"] = seqname
        by_accession[row["accession"]] = seqname
    return rows, columns, by_accession


def rewrite_fasta(input_path, output_path, seqnames):
    seen = set()
    count = 0
    with open(input_path) as source, open(output_path, "w") as output:
        for line in source:
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                accession = accession_from_fasta_id(identifier)
                if accession not in seqnames:
                    raise ValueError(
                        f"FASTA accession absent from metadata: {accession}"
                    )
                output.write(f">{seqnames[accession]}\n")
                seen.add(accession)
                count += 1
            else:
                output.write(line)

    missing = set(seqnames) - seen
    if missing:
        examples = ",".join(sorted(missing)[:10])
        raise ValueError(f"Metadata accessions absent from FASTA: {examples}")
    return count


def write_metadata(path, rows, columns):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def replace_files(fasta_path, metadata_path):
    rows, columns, seqnames = read_metadata(metadata_path)
    fasta_path = Path(fasta_path)
    metadata_path = Path(metadata_path)

    fasta_fd, fasta_temp = tempfile.mkstemp(
        prefix=f".{fasta_path.name}.",
        dir=fasta_path.parent,
    )
    metadata_fd, metadata_temp = tempfile.mkstemp(
        prefix=f".{metadata_path.name}.",
        dir=metadata_path.parent,
    )
    os.close(fasta_fd)
    os.close(metadata_fd)

    try:
        fasta_count = rewrite_fasta(fasta_path, fasta_temp, seqnames)
        write_metadata(metadata_temp, rows, columns)
        os.replace(fasta_temp, fasta_path)
        os.replace(metadata_temp, metadata_path)
    finally:
        for temp_path in (fasta_temp, metadata_temp):
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    return fasta_count, len(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()

    fasta_count, metadata_count = replace_files(args.fasta, args.metadata)
    print(f"Renamed FASTA records: {fasta_count}")
    print(f"Updated metadata rows: {metadata_count}")


if __name__ == "__main__":
    main()
