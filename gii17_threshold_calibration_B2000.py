#!/usr/bin/env python3
"""Calibrate distinguishable_rate threshold for GII.17 nextstrain clades.

Input: pre-aligned GII.17 genomes + clade assignment TSV (nextstrain-style).
Run A: clade-level majority-consensus refs (keep side).
Run B: per-strain refs (merge side, same-clade indistinguishability).
Same window/margin semantics as build_mapping_reference_group.py.
"""
import csv, re, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BT2 = "/Users/LuJ/Yunpan/software/Bioinf_software/bowtie2-2.3.3/bowtie2"
ALN = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/Phy_nextstrain_nf/data/GII17_genome/GII17.genome.alignment.fasta"
TSV = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/Phy_nextstrain_nf/data/GII17_genome/GII17.genome.assigned.tsv"
WORK = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/v2_test/gii17_threshold_calibration_B2000")
WINDOW, STEP, MIN_ALIGNED = 150, 25, 135
MIN_SCORE_MARGIN, MIN_NM_MARGIN = 6, 2
TARGET_MAX_MISMATCH = 0.10
BASES = frozenset("ACGT")

def ungapped(s): return s.replace("-", "")

def load_and_trim():
    records = list(SeqIO.parse(ALN, "fasta"))
    by_id = {r.id: str(r.seq).upper() for r in records}
    clade_of = {}
    dates = {}
    with open(TSV) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s, c = row["strain"].strip(), (row["clade"] or "").strip()
            if s in by_id and c:
                clade_of[s] = c
                dates[s] = (row["date"] or "")[:4]
    # User filter: clade B keeps only sequences collected >= 2000
    drop = {s for s, c in clade_of.items()
            if c == "B" and not (dates[s].isdigit() and int(dates[s]) >= 2000)}
    for s in drop:
        del clade_of[s]
    print(f"clade-B year filter: dropped {len(drop)} pre-2000 members")
    members = sorted(clade_of)
    arr = np.array([[c for c in by_id[m]] for m in members])
    keep = [j for j in range(arr.shape[1])
            if np.isin(arr[:, j], list("ACGT")).mean() >= 0.5]
    aligned = {m: "".join(by_id[m][j] for j in keep) for m in members}
    trimmed = [SeqRecord(Seq(aligned[m]), id=m, description="") for m in members]
    SeqIO.write(trimmed, WORK / "aligned_trimmed.fasta", "fasta")
    return aligned, clade_of, len(keep)

def consensus(ms, aligned):
    cols = list(zip(*[aligned[m] for m in ms]))
    out = []
    for col in cols:
        b = [x for x in col if x in BASES]
        out.append(Counter(b).most_common(1)[0][0] if b else "N")
    return "".join(out)

def make_windows(aligned, members):
    reads, meta, n = [], {}, 0
    for m in sorted(members):
        s = aligned[m]
        for st in range(0, len(s) - WINDOW + 1, STEP):
            frag = s[st:st + WINDOW]
            if any(c not in BASES for c in frag):
                continue
            n += 1
            rid = f"read_{n:08d}"
            meta[rid] = (m, st)
            reads.append(f">{rid}\n{frag}\n")
    (WORK / "windows.fasta").write_text("".join(reads))
    return meta

def run_bowtie2(refs, tag):
    d = WORK / tag
    d.mkdir(exist_ok=True)
    rf = d / "refs.fasta"
    rf.write_text("".join(f">{k}\n{v}\n" for k, v in refs.items()))
    idx = d / "index"
    subprocess.run([BT2 + "-build", str(rf), str(idx)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    sam = d / "win.sam"
    subprocess.run([BT2, "-x", str(idx), "-f", "-U", str(WORK / "windows.fasta"),
                    "-k", "2", "--very-sensitive-local", "--no-unal",
                    "--threads", "6", "-S", str(sam)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits = defaultdict(list)
    for line in sam.open():
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[5] == "*":
            continue
        aln = sum(int(l) for l, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5])
                  if op in "MI=X")
        tags = {}
        for t in f[11:]:
            p = t.split(":", 2)
            if len(p) == 3:
                tags[p[0]] = p[2]
        hits[f[0]].append({"ref": f[2], "as": int(tags.get("AS", -999)),
                           "nm": int(tags.get("NM", 999)), "len": aln})
    return hits

def evaluate(units, aligned, tag):
    refs = {u: ungapped(x["ref_aligned"]) for u, x in units.items()}
    assign = {m: u for u, x in units.items() for m in x["members"]}
    meta = make_windows(aligned, list(assign))
    hits = run_bowtie2(refs, tag)
    stats = defaultdict(Counter)
    for rid, (m, st) in meta.items():
        target = assign[m]
        ranked = sorted(hits.get(rid, []), key=lambda h: (-h["as"], h["nm"], h["ref"]))
        best = ranked[0] if ranked else None
        second = next((h for h in ranked[1:] if h["ref"] != best["ref"]), None) if best else None
        sm = best["as"] - second["as"] if best and second else 999
        nm = second["nm"] - best["nm"] if best and second else 999
        own = bool(best and best["ref"] == target)
        mapped = bool(best and best["len"] >= MIN_ALIGNED)
        t_slice = units[target]["ref_aligned"][st:st + WINDOW]
        comp = [(a, b) for a, b in zip(aligned[m][st:st + WINDOW], t_slice)
                if a in BASES and b in BASES]
        tm = bool(len(comp) >= MIN_ALIGNED and
                  sum(a != b for a, b in comp) / len(comp) <= TARGET_MAX_MISMATCH)
        dist = bool(mapped and own and (sm >= MIN_SCORE_MARGIN or nm >= MIN_NM_MARGIN))
        c = stats[target]
        c["total"] += 1; c["mapped"] += int(mapped); c["target"] += int(tm)
        c["own"] += int(own); c["dist"] += int(dist)
    return [{"unit": u, "windows": c["total"], "mapping_rate": c["mapped"]/c["total"],
             "target_rate": c["target"]/c["total"], "own_best_rate": c["own"]/c["total"],
             "distinguishable_rate": c["dist"]/c["total"]}
            for u, c in sorted(stats.items())]

def main():
    WORK.mkdir(parents=True, exist_ok=True)
    aligned, clade_of, aln_len = load_and_trim()
    by_clade = defaultdict(list)
    for m, c in clade_of.items():
        by_clade[c].append(m)
    print(f"members: {len(clade_of)}, clades: {len(by_clade)}, trimmed length: {aln_len}")
    for c in sorted(by_clade, key=lambda c: -len(by_clade[c])):
        print(f"  {c:15} {len(by_clade[c])}")

    unitsA = {c: {"ref_aligned": consensus(ms, aligned), "members": ms}
              for c, ms in by_clade.items()}
    rowsA = evaluate(unitsA, aligned, "runA_clade_panel")
    unitsB = {m: {"ref_aligned": aligned[m], "members": [m]} for m in clade_of}
    rowsB = evaluate(unitsB, aligned, "runB_strain_panel")

    vids = sorted(unitsA)
    D = np.zeros((len(vids), len(vids)))
    for i in range(len(vids)):
        for j in range(i + 1, len(vids)):
            a, b = unitsA[vids[i]]["ref_aligned"], unitsA[vids[j]]["ref_aligned"]
            ok = [(x, y) for x, y in zip(a, b) if x in BASES and y in BASES]
            D[i, j] = D[j, i] = sum(x != y for x, y in ok) / len(ok)

    for rows, fn in ((rowsA, "runA_clade_rates.tsv"), (rowsB, "runB_strain_rates.tsv")):
        with open(WORK / fn, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print("\n===== Run A: clade-level panel (KEEP side) =====")
    print(f"{'clade':15} {'n':>4} {'windows':>7} {'map':>6} {'target':>7} {'ownbest':>7} {'DIST':>7}")
    for r in sorted(rowsA, key=lambda r: r["distinguishable_rate"]):
        print(f"{r['unit']:15} {len(by_clade[r['unit']]):>4} {r['windows']:>7} "
              f"{r['mapping_rate']:>6.3f} {r['target_rate']:>7.3f} "
              f"{r['own_best_rate']:>7.3f} {r['distinguishable_rate']:>7.3f}")
    keep_min = min(r["distinguishable_rate"] for r in rowsA if r["windows"] > 0)
    print(f"KEEP-side minimum distinguishable_rate = {keep_min:.3f}")

    print("\n===== Run B: per-strain panel (MERGE side), by clade =====")
    byc = defaultdict(list)
    for r in rowsB:
        byc[r["unit"].split("|")[0] if "|" in r["unit"] else clade_of.get(r["unit"], "?")].append(r["distinguishable_rate"])
    for c in sorted(by_clade, key=lambda c: -len(by_clade[c])):
        v = sorted(byc.get(c, []))
        if not v:
            print(f"  {c:15} (no windows)")
            continue
        q = np.quantile(v, [0.1, 0.5, 0.9])
        print(f"  {c:15} n={len(v):>3}  p10={q[0]:.3f}  median={q[1]:.3f}  p90={q[2]:.3f}  min={v[0]:.3f} max={v[-1]:.3f}")

    print("\n===== clade consensus p-distance matrix =====")
    print("        " + " ".join(f"{v[:9]:>9}" for v in vids))
    for i, v in enumerate(vids):
        print(f"{v[:7]:7} " + " ".join(f"{D[i,j]:>9.4f}" for j in range(len(vids))))

    print("\n===== CONCLUSION =====")
    print(f"keep-side floor = {keep_min:.3f}")
    for T in (0.30, 0.50, 0.60, 0.80):
        print(f"  T={T:.2f}: 保留全部 clade={'PASS' if keep_min >= T else 'FAIL'}")

if __name__ == "__main__":
    main()
