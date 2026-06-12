#!/usr/bin/env python3
"""Download long Norovirus GI/GII records with a documented collection year."""

import argparse
import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GROUPS = {
    "gi": {"taxid": "122928", "label": "GI"},
    "gii": {"taxid": "122929", "label": "GII"},
}


def request(url, params, retries=6):
    params = dict(params)
    params["tool"] = "sewage_noro_ncbi_download"
    params["email"] = os.environ.get("NCBI_EMAIL", "norovirus-pipeline@example.org")
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    full_url = f"{url}?{urllib.parse.urlencode(params)}"

    for attempt in range(retries):
        try:
            with urllib.request.urlopen(full_url, timeout=120) as response:
                return response.read().decode("utf-8")
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(min(60, 2 ** attempt))


def esearch(taxid, min_length):
    term = f"txid{taxid}[Organism:exp] AND {min_length + 1}:100000[SLEN]"
    text = request(
        f"{EUTILS}/esearch.fcgi",
        {
            "db": "nuccore",
            "term": term,
            "retmode": "json",
            "retmax": 0,
            "usehistory": "y",
        },
    )
    result = json.loads(text)["esearchresult"]
    return term, int(result["count"]), result["webenv"], result["querykey"]


def iter_genbank_batches(webenv, query_key, count, batch_size):
    for start in range(0, count, batch_size):
        text = request(
            f"{EUTILS}/efetch.fcgi",
            {
                "db": "nuccore",
                "query_key": query_key,
                "WebEnv": webenv,
                "rettype": "gbwithparts",
                "retmode": "text",
                "retstart": start,
                "retmax": min(batch_size, count - start),
            },
        )
        yield start, text
        time.sleep(0.11 if os.environ.get("NCBI_API_KEY") else 0.36)


def parse_qualifier(block, name):
    match = re.search(
        rf'^\s+/{re.escape(name)}=(?:"((?:[^"]|"")*)"|([^\n]+))',
        block,
        flags=re.MULTILINE,
    )
    if not match:
        return ""
    value = (match.group(1) or match.group(2) or "").replace('""', '"')
    return " ".join(value.split())


def parse_genbank(text):
    for block in text.split("\n//"):
        if not block.strip():
            continue

        accession_match = re.search(r"^ACCESSION\s+(\S+)", block, re.MULTILINE)
        version_match = re.search(r"^VERSION\s+(\S+)", block, re.MULTILINE)
        definition_match = re.search(
            r"^DEFINITION\s+(.+?)(?=^ACCESSION\s)",
            block,
            flags=re.MULTILINE | re.DOTALL,
        )
        origin_match = re.search(r"^ORIGIN\s+(.+)$", block, re.MULTILINE | re.DOTALL)
        locus_match = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", block, re.MULTILINE)

        if not accession_match or not origin_match:
            continue

        sequence = re.sub(r"[^A-Za-z]", "", origin_match.group(1)).upper()
        accession = accession_match.group(1)
        version = version_match.group(1) if version_match else accession
        definition = ""
        if definition_match:
            definition = " ".join(definition_match.group(1).split())

        yield {
            "accession": accession,
            "version": version,
            "length": int(locus_match.group(1)) if locus_match else len(sequence),
            "definition": definition,
            "organism": parse_qualifier(block, "organism"),
            "collection_date": parse_qualifier(block, "collection_date"),
            "country": parse_qualifier(block, "country"),
            "isolate": parse_qualifier(block, "isolate"),
            "strain": parse_qualifier(block, "strain"),
            "sequence": sequence,
        }


def collection_year(value):
    years = [int(year) for year in re.findall(r"(?<!\d)(?:19|20)\d{2}(?!\d)", value)]
    return min(years) if years else None


def fasta_accessions(path):
    accessions = set()
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            identifier = line[1:].split()[0]
            match = re.search(
                r"(?:^|_)([A-Z]{1,4}_?\d{5,9})(?:\.\d+)?$",
                identifier,
            )
            if match:
                accessions.add(match.group(1))
    return accessions


def write_fasta_record(handle, record, label, year):
    accession = record["accession"]
    description = record["definition"].replace("\n", " ")
    handle.write(
        f">{label}_{accession} collection_date={record['collection_date']} "
        f"year={year} {description}\n"
    )
    sequence = record["sequence"]
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def download_group(group, output_dir, min_length, min_year, batch_size):
    config = GROUPS[group]
    term, candidate_count, webenv, query_key = esearch(config["taxid"], min_length)
    fasta_path = output_dir / f"ncbi_norovirus_{group}_gt{min_length}_collected_{min_year}_onward.fasta"
    metadata_path = output_dir / f"ncbi_norovirus_{group}_gt{min_length}_metadata.tsv"
    excluded_path = output_dir / f"ncbi_norovirus_{group}_gt{min_length}_excluded.tsv"

    columns = [
        "accession",
        "version",
        "group",
        "length",
        "collection_date",
        "collection_year",
        "organism",
        "country",
        "isolate",
        "strain",
        "definition",
    ]
    excluded_columns = ["accession", "length", "collection_date", "reason"]
    kept = 0
    missing_date = 0
    before_year = 0
    fetched = 0

    with (
        open(fasta_path, "w") as fasta_handle,
        open(metadata_path, "w", newline="") as metadata_handle,
        open(excluded_path, "w", newline="") as excluded_handle,
    ):
        metadata_writer = csv.DictWriter(metadata_handle, fieldnames=columns, delimiter="\t")
        excluded_writer = csv.DictWriter(excluded_handle, fieldnames=excluded_columns, delimiter="\t")
        metadata_writer.writeheader()
        excluded_writer.writeheader()

        for start, text in iter_genbank_batches(webenv, query_key, candidate_count, batch_size):
            for record in parse_genbank(text):
                fetched += 1
                year = collection_year(record["collection_date"])
                reason = ""
                if year is None:
                    reason = "missing_or_unparseable_collection_date"
                    missing_date += 1
                elif year < min_year:
                    reason = f"collection_year_before_{min_year}"
                    before_year += 1

                if reason:
                    excluded_writer.writerow(
                        {
                            "accession": record["accession"],
                            "length": record["length"],
                            "collection_date": record["collection_date"],
                            "reason": reason,
                        }
                    )
                    continue

                write_fasta_record(fasta_handle, record, config["label"], year)
                metadata_writer.writerow(
                    {
                        key: (
                            config["label"]
                            if key == "group"
                            else year
                            if key == "collection_year"
                            else record[key]
                        )
                        for key in columns
                    }
                )
                kept += 1

            print(f"{config['label']}: fetched {min(start + batch_size, candidate_count)}/{candidate_count}")

    return {
        "group": config["label"],
        "taxid": config["taxid"],
        "query": term,
        "candidate_count": candidate_count,
        "fetched_count": fetched,
        "kept_count": kept,
        "excluded_missing_or_unparseable_date": missing_date,
        "excluded_before_min_year": before_year,
        "fasta": str(fasta_path),
        "metadata": str(metadata_path),
        "excluded": str(excluded_path),
    }


def filter_genbank_file(group, genbank_path, output_dir, min_length, min_year):
    config = GROUPS[group]
    text = Path(genbank_path).read_text()
    records = list(parse_genbank(text))
    fasta_path = output_dir / f"ncbi_norovirus_{group}_gt{min_length}_collected_{min_year}_onward.fasta"
    metadata_path = output_dir / f"ncbi_norovirus_{group}_gt{min_length}_metadata.tsv"
    excluded_path = output_dir / f"ncbi_norovirus_{group}_gt{min_length}_excluded.tsv"
    columns = [
        "accession", "version", "group", "length", "collection_date",
        "collection_year", "organism", "country", "isolate", "strain",
        "definition",
    ]
    kept = []
    excluded = []

    for record in records:
        year = collection_year(record["collection_date"])
        if year is None:
            excluded.append((record, "missing_or_unparseable_collection_date"))
        elif year < min_year:
            excluded.append((record, f"collection_year_before_{min_year}"))
        else:
            kept.append((record, year))

    with open(fasta_path, "w") as handle:
        for record, year in kept:
            write_fasta_record(handle, record, config["label"], year)

    with open(metadata_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for record, year in kept:
            writer.writerow({
                key: (
                    config["label"] if key == "group"
                    else year if key == "collection_year"
                    else record[key]
                )
                for key in columns
            })

    with open(excluded_path, "w", newline="") as handle:
        columns_excluded = ["accession", "length", "collection_date", "reason"]
        writer = csv.DictWriter(handle, fieldnames=columns_excluded, delimiter="\t")
        writer.writeheader()
        for record, reason in excluded:
            writer.writerow({
                "accession": record["accession"],
                "length": record["length"],
                "collection_date": record["collection_date"],
                "reason": reason,
            })

    return {
        "group": config["label"],
        "taxid": config["taxid"],
        "source_genbank": str(genbank_path),
        "candidate_count": len(records),
        "fetched_count": len(records),
        "kept_count": len(kept),
        "excluded_missing_or_unparseable_date": sum(
            reason == "missing_or_unparseable_collection_date"
            for _, reason in excluded
        ),
        "excluded_before_min_year": sum(
            reason.startswith("collection_year_before_")
            for _, reason in excluded
        ),
        "fasta": str(fasta_path),
        "metadata": str(metadata_path),
        "excluded": str(excluded_path),
    }


def compare_with_current(summary, current_path, output_dir):
    downloaded = fasta_accessions(summary["fasta"])
    current = fasta_accessions(current_path)
    group = summary["group"].lower()
    overlap = sorted(downloaded & current)
    new_accessions = sorted(downloaded - current)
    missing_accessions = sorted(current - downloaded)
    list_paths = {
        "overlap_accession_list": output_dir / f"{group}_overlap_accessions.txt",
        "new_accession_list": output_dir / f"{group}_new_in_ncbi_download.txt",
        "current_missing_accession_list": output_dir / f"{group}_current_not_in_filtered_download.txt",
    }
    for key, path in list_paths.items():
        values = (
            overlap if key == "overlap_accession_list"
            else new_accessions if key == "new_accession_list"
            else missing_accessions
        )
        path.write_text("\n".join(values) + ("\n" if values else ""))

    summary["current_dataset"] = str(current_path)
    summary["current_accessions"] = len(current)
    summary["overlap_accessions"] = len(overlap)
    summary["new_in_ncbi_download"] = len(new_accessions)
    summary["current_not_in_filtered_download"] = len(missing_accessions)
    summary["current_coverage_percent"] = round(
        100 * len(overlap) / len(current), 2
    ) if current else 0.0
    summary["download_overlap_percent"] = round(
        100 * len(overlap) / len(downloaded), 2
    ) if downloaded else 0.0
    summary.update({key: str(path) for key, path in list_paths.items()})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="ncbi_download_2026-06-11")
    parser.add_argument("--min-length", type=int, default=5800)
    parser.add_argument("--min-year", type=int, default=2000)
    parser.add_argument("--batch-size", type=int, default=200)
    parser.add_argument("--gi-current", default="gi_with_genotype-original.fasta")
    parser.add_argument("--gii-current", default="gii_with_genotype-original.fasta")
    parser.add_argument("--gi-genbank", help="Use an already downloaded GI GenBank file")
    parser.add_argument("--gii-genbank", help="Use an already downloaded GII GenBank file")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for group, current_path in (("gi", args.gi_current), ("gii", args.gii_current)):
        genbank_path = args.gi_genbank if group == "gi" else args.gii_genbank
        if genbank_path:
            summary = filter_genbank_file(
                group, genbank_path, output_dir, args.min_length, args.min_year
            )
        else:
            summary = download_group(
                group, output_dir, args.min_length, args.min_year, args.batch_size
        )
        if Path(current_path).exists():
            compare_with_current(summary, current_path, output_dir)
        summaries.append(summary)

    summary_path = output_dir / "download_and_comparison_summary.json"
    summary_path.write_text(json.dumps(summaries, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summaries, indent=2, ensure_ascii=False))
    print(f"Summary written to {summary_path}")


if __name__ == "__main__":
    main()
