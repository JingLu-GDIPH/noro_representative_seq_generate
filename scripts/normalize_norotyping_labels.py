#!/usr/bin/env python3
"""Normalize genotype labels missed by the local norotyping header parser."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


UNKNOWN = {"", "Unknown", "NA", "N/A", "Unassigned"}


def genotype_from_reference(reference: str) -> str | None:
    parts = (reference or "").split()
    if not parts:
        return None
    token = parts[0]
    match = re.search(
        r"(?:^|[_|])(G?(?:IX|VIII))\.(NA\d+|\d+)(?:[_|.]|$)",
        token,
        re.IGNORECASE,
    )
    if match:
        group = match.group(1).upper()
        if not group.startswith("G"):
            group = f"G{group}"
        return f"{group}.{match.group(2).upper()}"
    return None


def normalize_row(row: dict[str, str]) -> bool:
    changed = False

    vp1 = row.get("vp1_genotype", "")
    inferred_vp1 = genotype_from_reference(row.get("vp1_nearest_reference", ""))
    if inferred_vp1 and vp1 in UNKNOWN:
        row["vp1_genotype"] = inferred_vp1
        changed = True

    if "name_genotype" in row and "vp1_matches_name" in row:
        name = row.get("name_genotype", "").strip()
        vp1 = row.get("vp1_genotype", "").strip()
        row["vp1_matches_name"] = str(bool(name and vp1 and name == vp1)).lower()

    if "rdrp_vp1_combo" in row:
        rdrp = row.get("rdrp_genotype", "").strip() or "Unknown"
        vp1 = row.get("vp1_genotype", "").strip() or "Unknown"
        row["rdrp_vp1_combo"] = f"{rdrp}/{vp1}"

    return changed


def normalize_table(path: Path) -> int:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    changed_count = sum(1 for row in rows if normalize_row(row))
    if changed_count:
        backup = Path(str(path) + ".before_normalize_gix_labels")
        if not backup.exists():
            backup.write_text(path.read_text())
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    return changed_count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", nargs="+", type=Path)
    args = parser.parse_args()

    for table in args.tables:
        changed = normalize_table(table)
        print(f"{table}: normalized {changed} rows")


if __name__ == "__main__":
    main()
