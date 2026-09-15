#!/usr/bin/env python3
"""Audit final mapping references against the pipeline's existing QC rules."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from Bio import AlignIO, SeqIO


BASES = frozenset("ACGT")


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def majority_sequence(member_ids: list[str], aligned: dict[str, str]) -> str:
    output = []
    for column in zip(*(aligned[member_id] for member_id in member_ids)):
        counts = Counter(base for base in column if base in BASES)
        output.append(counts.most_common(1)[0][0] if counts else "N")
    return "".join(output)


def private_positions(reference: str, query: str) -> list[int]:
    return [
        position
        for position, (left, right) in enumerate(zip(reference, query), start=1)
        if left in BASES and right in BASES and left != right
    ]


def count_snp_clusters(
    positions: list[int], window: int = 101, cutoff: int = 6
) -> int:
    """Use the same SNP-cluster implementation as the mapping-reference branch."""
    if not positions:
        return 0
    left = 0
    in_cluster = False
    clusters = 0
    for right, position in enumerate(positions):
        while positions[left] < position - window + 1:
            left += 1
        if right - left + 1 > cutoff:
            if not in_cluster:
                clusters += 1
                in_cluster = True
        else:
            in_cluster = False
    return clusters


def longest_run(sequence: str, characters: str) -> int:
    pattern = f"[{re.escape(characters)}]+"
    return max((len(match.group()) for match in re.finditer(pattern, sequence)), default=0)


def audit_group(result_dir: Path, group: str) -> list[dict]:
    group_dir = result_dir / group
    final_fasta = result_dir / f"{group}_mapping_references.fasta"
    final_records = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(final_fasta, "fasta")
    }
    rows = []

    for reference_dir in sorted((group_dir / "06_group_references").glob("mapping_ref_*")):
        group_key = reference_dir.name.removeprefix("mapping_ref_")
        alignment_path = group_dir / "05_iqtree" / f"aligned_{group_key}.fasta"
        alignment = AlignIO.read(alignment_path, "fasta")
        aligned = {
            record.id: str(record.seq).upper()
            for record in alignment
        }
        alignment_length = alignment.get_alignment_length()
        group_root = majority_sequence(sorted(aligned), aligned)

        members = defaultdict(list)
        for row in read_tsv(reference_dir / "final_membership.tsv"):
            members[row["reference"]].append(row["raw_sequence"])

        raw_qc = {
            row["sequence_id"]: row
            for row in read_tsv(reference_dir / "raw_reference_qc.tsv")
        }
        validation = {
            row["reference"]: row
            for row in read_tsv(reference_dir / "final_window_validation.tsv")
        }

        for source in read_tsv(reference_dir / "final_reference_sources.tsv"):
            output_id = source["output_id"]
            internal_id = source["internal_reference"]
            final_sequence = final_records[output_id]
            member_ids = members[internal_id]

            if source["source_type"] in {"observed_raw", "medoid"}:
                aligned_reference = aligned[source["source_record"]]
                if aligned_reference.replace("-", "") != final_sequence:
                    raise ValueError(f"Source sequence mismatch for {output_id}")
            elif len(final_sequence) == alignment_length:
                aligned_reference = final_sequence
            else:
                raise ValueError(
                    f"Cannot restore alignment columns for {output_id}: "
                    f"{len(final_sequence)} != {alignment_length}"
                )

            local_root = majority_sequence(sorted(member_ids), aligned)
            local_positions = private_positions(local_root, aligned_reference)
            group_positions = private_positions(group_root, aligned_reference)

            source_qc = raw_qc.get(source["source_record"], {})
            source_flag = source_qc.get("flag", "") if source["source_record"] else ""
            source_status = source_flag or (
                "PASS" if source["source_record"] else "NOT_APPLICABLE"
            )

            flagged_members = sum(
                raw_qc.get(member_id, {}).get("flag") == "SNP_CLUSTERS>1"
                for member_id in member_ids
            )
            unassessable_members = sum(
                raw_qc.get(member_id, {}).get("flag")
                in {"SINGLETON_NO_REFERENCE", "SMALL_GROUP_NO_STABLE_REFERENCE"}
                for member_id in member_ids
            )

            other_codes = sorted(set(final_sequence) - set("ACGTN-"))
            n_gap_percent = (
                100.0
                * (final_sequence.count("N") + final_sequence.count("-"))
                / len(final_sequence)
            )
            validation_row = validation[internal_id]

            notes = []
            status = "PASS"
            if other_codes:
                status = "FAIL"
                notes.append("non_ACGTN_base")
            if n_gap_percent > 10.0:
                status = "FAIL"
                notes.append("N_gap_gt_10pct")
            if source_flag == "SNP_CLUSTERS>1":
                status = "FAIL"
                notes.append("flagged_raw_source")

            local_clusters = (
                count_snp_clusters(local_positions) if len(member_ids) > 1 else None
            )
            if status == "PASS" and local_clusters is not None and local_clusters > 1:
                status = "REVIEW"
                notes.append("local_SNP_clusters_gt_1")
            if (
                status == "PASS"
                and source_flag
                in {"SINGLETON_NO_REFERENCE", "SMALL_GROUP_NO_STABLE_REFERENCE"}
            ):
                status = "UNASSESSABLE"
                notes.append("no_stable_raw_QC_reference")

            rows.append(
                {
                    "group": group.upper(),
                    "reference_id": output_id,
                    "source_type": source["source_type"],
                    "source_record": source["source_record"],
                    "group_size": len(aligned),
                    "assigned_members": len(member_ids),
                    "length": len(final_sequence),
                    "n_gap_percent": f"{n_gap_percent:.6f}",
                    "other_ambiguity_count": sum(
                        base not in set("ACGTN-") for base in final_sequence
                    ),
                    "other_ambiguity_codes": ",".join(other_codes),
                    "longest_n_run": longest_run(final_sequence, "N"),
                    "raw_source_qc": source_status,
                    "flagged_input_members": flagged_members,
                    "unassessable_input_members": unassessable_members,
                    "local_private_snps": (
                        len(local_positions) if len(member_ids) > 1 else ""
                    ),
                    "local_snp_clusters": (
                        local_clusters if local_clusters is not None else ""
                    ),
                    "group_private_snps": len(group_positions),
                    "group_snp_clusters": count_snp_clusters(group_positions),
                    "target_mapping_rate": validation_row["target_mapping_rate"],
                    "distinguishable_rate": validation_row["distinguishable_rate"],
                    "mapping_decision": validation_row["decision"],
                    "audit_status": status,
                    "audit_notes": ";".join(notes),
                }
            )

    missing = set(final_records) - {row["reference_id"] for row in rows}
    if missing:
        raise ValueError(f"References missing from audit: {sorted(missing)}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    for group in ("gi", "gii"):
        rows.extend(audit_group(args.result_dir, group))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    for group in ("GI", "GII"):
        group_rows = [row for row in rows if row["group"] == group]
        statuses = Counter(row["audit_status"] for row in group_rows)
        print(group, len(group_rows), dict(statuses))


if __name__ == "__main__":
    main()
