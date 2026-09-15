#!/usr/bin/env python3
"""Parse GenBank reference records and append them to the reference FASTA +
gene-coordinate TSV so the gene-annotation pipeline can use them.

Two norovirus genotypes (GII.9, GIX.1) were missing from the reference set.
Their GenBank records (``reference/gii9.gb``, ``reference/gix.gb``) provide both
the sequence and the ORF1 mat_peptide / ORF2 / ORF3 feature coordinates.

For each GenBank record this script emits:
  * a FASTA record ``<VP1genotype>_<accession>`` appended to the relevant
    group's reference FASTA (gi_with_genotype.fasta / gii_with_genotype.fasta);
  * rows appended to ``reference/genome_gene_coordinates.tsv`` mirroring the
    existing schema, translating GenBank feature qualifiers to the TSV columns:

      protein  <- feature product, mapped to the TSV's canonical names:
                     p48/NTPase/p22/VPg kept; "Pro" -> "3CLpro";
                     ORF1 from the parent CDS; VP1/VP2 from ORF2/ORF3 CDS.
      start/end <- feature location (1-based inclusive, join-aware)
      strand   <- "+" (norovirus genomes are plus strand)
      source   <- "genbank" for CDS-level (ORF1/VP1/VP2),
                  "orf1_relative" for mat_peptide sub-units.

Idempotent: re-running drops any previously-appended rows for the GenBank
accessions before re-adding them.
"""

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO


#: Map GenBank mat_peptide / CDS product labels to TSV canonical protein names.
PRODUCT_MAP = {
    "p48": "p48",
    "NTPase": "NTPase",
    "p22": "p22",
    "VPg": "VPg",
    "Pro": "3CLpro",      # 3C-like protease
    "3CLpro": "3CLpro",
    "RdRp": "RdRp",
    "VP1": "VP1",
    "VP2": "VP2",
}

#: VP1 genotype label per accession, derived from the GenBank /note.
ACCESSION_TO_VP1 = {
    "MW305520": ("GII", "GII.9"),
    "MN227775": ("GII", "GIX.1"),
}


def location_bounds(feature):
    """Return (start, end) 1-based inclusive, handling simple join locations."""
    parts = feature.location.parts
    start = min(int(p.start) + 1 for p in parts)
    end = max(int(p.end) for p in parts)
    return start, end


def feature_rows(record, group, vp1):
    """Yield gene-coordinate dicts for one GenBank record."""
    accession = record.id.split(".")[0]
    seq_len = len(record.seq)
    seen_orf1_cds = False
    for feature in record.features:
        ftype = feature.type
        product = feature.qualifiers.get("product", [""])[0]
        gene = feature.qualifiers.get("gene", [""])[0]

        # ORF1 as a whole (from the parent CDS).
        if ftype == "CDS" and gene == "ORF1":
            start, end = location_bounds(feature)
            yield {
                "group": group, "accession": accession, "genotype": vp1,
                "sequence_length": seq_len, "protein": "ORF1",
                "start": start, "end": end, "strand": "+", "source": "genbank",
            }
            seen_orf1_cds = True
            continue

        # ORF2 -> VP1, ORF3 -> VP2 (major / minor capsid).
        if ftype == "CDS" and gene in ("ORF2", "ORF3"):
            protein = PRODUCT_MAP.get(product, product)
            start, end = location_bounds(feature)
            yield {
                "group": group, "accession": accession, "genotype": vp1,
                "sequence_length": seq_len, "protein": protein,
                "start": start, "end": end, "strand": "+", "source": "genbank",
            }
            continue

        # ORF1 mature peptides (sub-units).
        if ftype == "mat_peptide":
            protein = PRODUCT_MAP.get(product, product)
            start, end = location_bounds(feature)
            yield {
                "group": group, "accession": accession, "genotype": vp1,
                "sequence_length": seq_len, "protein": protein,
                "start": start, "end": end, "strand": "+",
                "source": "orf1_relative",
            }


def integrate(genbank_files, reference_fastas, coords_tsv):
    """Append GenBank-derived FASTA + coord rows; idempotent per accession."""
    coords_path = Path(coords_tsv)
    new_accessions = set()
    new_rows = []
    fasta_appends = {g: [] for g in ("GI", "GII")}

    for gb_path in genbank_files:
        for record in SeqIO.parse(gb_path, "genbank"):
            accession = record.id.split(".")[0]
            if accession not in ACCESSION_TO_VP1:
                print(f"[WARN] {gb_path}: accession {accession} not in "
                      f"ACCESSION_TO_VP1 map; skipping")
                continue
            group, vp1 = ACCESSION_TO_VP1[accession]
            new_accessions.add(accession)
            # FASTA record: <VP1genotype>_<accession>
            fasta_id = f"{vp1}_{accession}"
            fasta_appends[group].append((fasta_id, str(record.seq)))
            new_rows.extend(feature_rows(record, group, vp1))
            print(f"{gb_path}: {accession} group={group} vp1={vp1} "
                  f"len={len(record.seq)} features={sum(1 for _ in feature_rows(record, group, vp1))}")

    # Rewrite coords TSV: drop any prior rows for these accessions, then append.
    with open(coords_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        kept = [row for row in reader if row["accession"] not in new_accessions]
    dropped = sum(1 for r in kept if r["accession"] in new_accessions)
    print(f"\nCoords TSV: {len(kept)} kept (dropped {dropped} prior rows for "
          f"re-added accessions), appending {len(new_rows)} new rows")

    with open(coords_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)
        for row in new_rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    # Append FASTA records (idempotent: drop existing lines with these IDs first).
    for group, path in reference_fastas.items():
        if not fasta_appends.get(group):
            continue
        existing = []
        new_ids = {fid for fid, _ in fasta_appends[group]}
        with open(path) as handle:
            block = []
            current_id = None
            for line in handle:
                if line.startswith(">"):
                    if current_id is not None:
                        existing.append((current_id, block))
                    current_id = line[1:].split()[0]
                    block = [line]
                else:
                    block.append(line)
            if current_id is not None:
                existing.append((current_id, block))
        existing = [(cid, blk) for cid, blk in existing if cid not in new_ids]
        with open(path, "w") as handle:
            for _cid, blk in existing:
                handle.writelines(blk)
            for fid, seq in fasta_appends[group]:
                handle.write(f">{fid}\n")
                for i in range(0, len(seq), 60):
                    handle.write(seq[i:i + 60] + "\n")
        print(f"FASTA {path}: now {len(existing) + len(fasta_appends[group])} records "
              f"(appended {len(fasta_appends[group])})")


def main():
    parser = argparse.ArgumentParser(
        description="Integrate GenBank reference records into the reference "
                    "FASTA + gene-coordinate TSV used by the annotation pipeline.")
    parser.add_argument("--genbank", nargs="+", required=True,
                        help="GenBank reference file(s), e.g. reference/gii9.gb reference/gix.gb")
    parser.add_argument("--reference-fastas", nargs="+", required=True,
                        help="group=path pairs, e.g. GII=reference/gii_with_genotype.fasta")
    parser.add_argument("--coords", required=True,
                        help="reference/genome_gene_coordinates.tsv")
    args = parser.parse_args()

    reference_fastas = {}
    for item in args.reference_fastas:
        if "=" not in item:
            sys.exit(f"--reference-fastas expects GROUP=PATH, got {item}")
        group, path = item.split("=", 1)
        reference_fastas[group] = path

    integrate(args.genbank, reference_fastas, args.coords)


if __name__ == "__main__":
    main()
