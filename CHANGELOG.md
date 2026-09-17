# Changelog

## v2.1.1 (2026-09-15)

- Renamed `main.nf` -> `main_probe.nf` for clarity (probe-consensus pipeline);
  run with `nextflow run main_probe.nf ...`. Docs and config comments updated.


## v2.1.0 (2026-09-15)

### Repository cleanup — two production pipelines only
- Scripts reorganized into `scripts/common/` (shared), `scripts/probe_consensus/`
  (pipeline 1), and `scripts/mapping_reference/` (pipeline 2); `main.nf` and
  `mapping_reference.nf` path references updated.
- Removed data-preparation chain (NCBI download/norotyping/renaming), one-off
  calibration & validation scripts, plotting/QC extras, `tests/`, and stale docs
  (METHODS_AND_PIPELINE.md, MAPPING_REFERENCE_README.md, SCRIPT_VERSIONS.md).
  All recoverable from git history (v2.0.0 commit).
- README rewritten for the cleaned layout; STRAIN_PANEL_METHOD.md paths updated.
- `result_paper/` untracked (kept locally; embedded in the analysis docx).


## v2.0.0 (2026-09-15)

### Strain-level medoid panel (final recipe)
- `--representative medoid`: references are now REAL observed sequences (member
  closest by patristic distance to the clade root node), with QC-preferred
  selection and `MEDOIDQC` fallback tagging. Hybrid ASR path kept behind
  `--representative asr`.
- Strain-level selection: `--max-within-d95 0.005 --min-parent-d95 0.005
  --min-parent-delta 0.001 --max-within-max 0.02`; UFBoot gate dropped
  (`--min-ufboot 0`) — SH-aLRT-only clades.
- VP1 segment-level merge veto (`--veto-block-windows 15 --veto-window-snv 2
  --veto-min-window-coverage 0.5`): a divergent block of >=15 consecutive VP1
  windows (>=2 SNVs each) forbids merging; VP1 located by annotation transfer
  (`mafft --add --keeplength`); partial-genome coverage guard included.
- Batch merge of disjoint pairs (replaces one-pair-per-iteration) and batch
  split; MRCA cascade merge removed (240->2 avalanche root cause).
- Umbrella split routing (own_best < 0.5) with `--umbrella-gate structural`
  (>=2 child clades each holding >= min-tips of the unit's own members) or
  `members` (legacy >=25).
- Oscillation detection (member-partition signature) + best-state snapshot
  (score: fewest flagged units, then fewer references) emitted as final output.
- Final shadow sweep (`--shadow-sweep-distance 0.005`): collapse residual
  <=0.5%-apart reference pairs after the loop, revalidated in `final_sweep/`.
- Fixed: split no longer pulls members owned by other units; restore/sweep now
  runs before output writing (output consistency).

### Validation
- GII.P16_GII.4 (504 seqs): 37 refs, all real, min spacing 0.64%, zero shadows.
- GII.17 (7 groups, 664-seq main): 34 refs; Kawasaki_308/Kawasaki_323/Romania
  each receive dedicated medoid references (Node276 collapse resolved).
- Forward derivation (derive_params_forward.py): R99=0.805 for 99%-similar
  pairs; mixed-strain recovery 0.000% base error / 100% discriminating sites /
  0% chimera down to 10% minority with a distinguished panel.


## v1.4.0 (2026-07-01)

### Added
- `ORIENT_SEQUENCES` process and `scripts/normalize_orientation.py`: before MAFFT,
  every cluster is put on a single strand. Each record is compared (by 15-mer
  overlap, forward vs reverse complement) against the cluster's longest A/C/G/T
  record, and any record whose reverse complement matches the reference better is
  reverse-complemented in place. A full audit of the GI/GII inputs found three
  reverse-strand records in two GII clusters (`GII.P31_GII.4_KX158285`,
  `GII.P7_GII.6_KX158282`, `GII.P7_GII.6_KX158286`); without this step the
  reversed records were forced into MAFFT on the wrong strand, producing ~1,500
  spurious alignment gaps and meaningless pairwise similarities.
- `tests/test_normalize_orientation.py`: unit tests covering reverse-strand
  flipping, forward-record preservation, single-record passthrough, N/gap-padded
  input handling, and the minimum-overlap guard.
- `scripts/analyze_raw_representative_qc.py`: nextclade-style QC of raw-sequence
  representatives against their genotype clade-root consensus. Counts
  `private_mutations` and `snp_clusters` (101 nt window, >6 private SNPs) and
  tags suspects `SNP_CLUSTERS>1` / `HIGH_PRIVATE`. Singleton clusters are tagged
  `SINGLETON_NO_REFERENCE` (undecidable, not a QC pass).
- QC-clean representative sets `gi/gii_all_final_consensus_filtered.fa` in
  `07_final_results/`, with raw reps tagged `SNP_CLUSTERS>1` removed
  (GI 42 -> 29, GII 117 -> 92).
- `results_new/qc_analysis/REPORT.md`: method, results, and rep-vs-sibling
  baseline cross-check.

### Changed
- `main.nf` workflow now chains `SPLIT_CLUSTERS -> ORIENT_SEQUENCES -> ALIGN`.
- `CONSENSUS` filename parsing no longer anchors on the leading `cluster__`, so
  it still extracts genotype/cluster id from the `oriented_`-prefixed filenames.
- `nextflow.config` manifest version bumped to `1.4.0` and adds lightweight
  resources for the `ORIENT_SEQUENCES` process.
- `README.md`, `SCRIPT_VERSIONS.md`: document the representative QC step, the
  `_filtered.fa` outputs, and recommend `_filtered.fa` for downstream use.

## v1.3.1 (2026-06-22)

### Changed
- MAFFT post-processing now removes every alignment column with A/C/G/T
  coverage below 50%; `N` and gap characters are not counted as valid coverage.
- Final representative FASTA records are written without alignment gaps, so
  low-coverage internal insertion columns cannot inflate consensus sequence
  length.
- Consensus base calling now counts only A/C/G/T bases; `N` and gaps are
  excluded from majority consensus decisions.

### Results
- Refreshed formal outputs in `results_new/gi` and `results_new/gii` using the
  norotyping-renamed GI/GII FASTA inputs.
- GII final representatives now have no gap characters and no sequences above
  8 kb; GII.4 length range is 6038-7538 bp and GII.6 length range is
  5030-7545 bp.

## v1.2.1 (2026-06-15)

### Added
- `analyze_consensus_representativeness.py`: adds VP1/RdRp labels to final
  consensus FASTA IDs, splits downloaded and consensus sequences by VP1 genotype,
  runs MAFFT/FastTree per genotype, and reports representativeness metrics.
- `plot_consensus_representativeness_trees.py`: renders selected per-genotype
  raw+consensus trees with consensus tips highlighted.

### Changed
- Default runtime thread budget is now `--threads 8`, matching the recommended
  setting for this 10-logical-core workstation.
- `nextflow.config` now derives local executor CPUs, queue size, per-task CPUs,
  and `maxForks` from `--threads` to improve CPU use while avoiding accidental
  oversubscription.
- `CLUSTER` and final global pruning can use the full thread budget; MAFFT and
  IQ-TREE based tasks use up to 4 threads per task with controlled concurrency.
- `main.nf` validates that `--threads` is a positive integer.

## v1.2.0 (2026-06-15)

### Added
- `filter_vcluster_pairs_by_genotype.py`: filters vclust ANI pairs so clustering
  is performed only within the same VP1 genotype (`GI.x`, `GII.x`, or `GIX.x`).
- `README.md` and `SCRIPT_VERSIONS.md` for reproducible use by other scripts
  and downstream workflows.

### Changed
- Default clustering parameters now use ANI `0.90`, qcov `0.65`, rcov `0.65`,
  len_ratio `0.65`, and complete linkage.
- MAFFT alignment is always run for multi-sequence clusters after removing
  pre-existing gaps, so padded FASTA records are not mistaken for aligned data.
- Terminal alignment columns are trimmed when A/C/G/T coverage is below 50%;
  `N` is not counted as valid coverage.
- Main workflow labels and configuration comments now refer to `main.nf` and
  support both GI and GII datasets.

### Results
- Formal outputs are kept in `results_new/gi` and `results_new/gii`.
- The refreshed GII run produced 87 VP1 genotype-pure clusters, 24 singleton
  clusters, 128 final representative sequences, and a maximum final similarity
  of 94.87%.

## v1.1.0 (2026-06-12)

### Fixed
- `prune_redundant_consensus.py`: added robust error handling for MAFFT failures
  (empty output, non-zero exit, dropped sequence IDs) instead of silently
  crashing with a bare `KeyError`.
- `COLLECT_RESULTS` step now correctly calls `prune_redundant_consensus.py` to
  remove cross-cluster redundant sequences before final output. Previous runs
  using cached Nextflow work dirs may still use the old command; clear the
  `work/` directory or run without `-resume` to pick up the fix.

### Added
- `COLLECT_RESULTS` process: cross-cluster redundancy pruning via
  `prune_redundant_consensus.py` (writes intermediate
  `all_final_consensus.raw.fasta`, then prunes to `all_final_consensus.fasta`).
- `.gitignore`, `VERSION`, `CHANGELOG.md` for version control.

### Changed
- `nextflow.config` manifest: name `noro-consensus`, version `1.1.0`.

---

## v1.0.0 (2026-06-11)

### Initial release
- Eight-step Nextflow pipeline: quality filter → vclust clustering → split →
  MAFFT alignment → phylogenetic consensus (IQ-TREE ASR) → enhanced iterative
  validation → collect results → final validation.
- Supports GI and GII norovirus datasets.
- Python scripts for consensus generation, trimming, validation, and plotting.
- `METHODS_AND_PIPELINE.md`: full method description (Chinese + English).
