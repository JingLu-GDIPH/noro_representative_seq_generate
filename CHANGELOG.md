# Changelog

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
