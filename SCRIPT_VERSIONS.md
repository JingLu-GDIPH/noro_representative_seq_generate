# Script Version Notes

Pipeline version: `1.4.0`

Last updated: 2026-07-01

## Workflow Entry Points

| File | Role | Version notes |
| --- | --- | --- |
| `main.nf` | Main Nextflow workflow | v1.4.0 inserts an `ORIENT_SEQUENCES` step between cluster splitting and MAFFT so that records stored on the reverse strand (e.g. `GII.P31_GII.4_KX158285`) are reverse-complemented before alignment; `CONSENSUS` filename parsing now tolerates the `oriented_` prefix. v1.3.1 defaults to VP1 genotype-stratified grouping, MAFFT realignment, full-column A/C/G/T coverage masking at 0.50, iterative consensus validation, final gap-free representative output, and positive `--threads` validation. |
| `nextflow.config` | Local runtime config | v1.4.0 manifest version and `ORIENT_SEQUENCES` process resources. Local executor, conda environment, reports, traces, and dynamic process resources derived from `--threads`. |
| `environment.yml` | Conda environment | Defines the `noro-consensus` software environment used by workflow processes. |

## Data Preparation Scripts

| Script | Role | Called by workflow |
| --- | --- | --- |
| `scripts/download_ncbi_norovirus_by_collection_date.py` | Download NCBI GI/GII nucleotide records with sequence length >5800 bp and documented collection date after 2000. | No |
| `scripts/type_vp1_and_update_metadata.py` | Assign VP1 genotype using the IPHnano norovirus genotyping references and update metadata. | No |
| `scripts/rename_vp1_typed_sequences.py` | Rename FASTA records to `genotype_accession_collectiondate` and add `seqname` to metadata. | No |
| `scripts/run_consensus_genotyping.py` | Run the IPHnano dual-region genotyper on final GI/GII consensus FASTA files with multi-worker execution and a MAFFT `--quiet` wrapper. | No |
| `scripts/analyze_consensus_representativeness.py` | Add VP1/RdRp labels to consensus FASTA IDs, split raw and consensus sequences by VP1 genotype, run MAFFT/FastTree per genotype, and calculate representativeness metrics. | No |
| `scripts/plot_consensus_representativeness_trees.py` | Render selected raw+consensus genotype trees with consensus tips highlighted. | No |

## Main Pipeline Scripts

| Script | Role | Called by workflow |
| --- | --- | --- |
| `scripts/filter_n_sequences.py` | Remove records where `N + gap` percentage exceeds the quality threshold. | Yes, `MERGE_FASTA` |
| `scripts/filter_vcluster_pairs_by_genotype.py` | Keep only vclust ANI pairs where query and reference have the same VP1 genotype. | Yes, `CLUSTER` |
| `scripts/process_vcluster_results.py` | Convert vclust clusters into `cluster_info.csv` and clustered FASTA. | Yes, `CLUSTER` |
| `scripts/split_clusters.py` | Split clustered sequences into one FASTA per cluster. | Yes, `SPLIT_CLUSTERS` |
| `scripts/normalize_orientation.py` | Reverse-complement any record stored opposite to the cluster majority strand, using k-mer overlap against the longest A/C/G/T record as the orientation reference, so MAFFT receives a consistently oriented cluster. | Yes, `ORIENT_SEQUENCES` |
| `scripts/trim_alignment_ends.py` | Remove all alignment columns with A/C/G/T coverage below threshold; `N` and gap are invalid coverage. | Yes, `ALIGN` |
| `scripts/generate_consensus.py` | Build tree/ASR and generate first-pass consensus sequences from phylogenetic boundary nodes. | Yes, `CONSENSUS` |
| `scripts/enhanced_consensus_with_validation.py` | Iteratively validate and refine per-cluster consensus output until internal similarity threshold is met or max iterations are reached. | Yes, `ENHANCED_CONSENSUS` |
| `scripts/trim_sequence_ends.py` | Remove terminal `N`/gap characters and, when requested by the workflow, all alignment gaps from final representative sequences. | Yes, `ENHANCED_CONSENSUS` |
| `scripts/combine_consensus_fastas.py` | Merge per-cluster FASTA files and keep identifiers globally unique. | Yes, `COLLECT_RESULTS` |
| `scripts/prune_redundant_consensus.py` | Globally remove final representatives with similarity above threshold after MAFFT alignment. | Yes, `COLLECT_RESULTS` |
| `scripts/validate_consensus_internal.py` | Validate final pairwise similarities and write the report. | Yes, `FINAL_VALIDATION` |

## Optional Validation and Plotting Scripts

| Script | Role |
| --- | --- |
| `scripts/analyze_raw_representative_qc.py` | Flag raw-sequence representatives with excessive private mutations / SNP clusters vs their genotype clade root (nextclade-style), to detect probable sequencing errors. Outputs a TSV used to build `*_all_final_consensus_filtered.fa`. |
| `scripts/validate_representativeness.py` | Exploratory representativeness validation. |
| `scripts/validate_representativeness_blast.py` | BLAST-based representativeness validation. |
| `scripts/validate_representativeness_final.py` | MAFFT-based cluster consensus representativeness validation. |
| `scripts/validate_representativeness_mafft.py` | MAFFT-based original-to-consensus validation. |
| `scripts/plot_validation_results.py` | Generate validation plots. |
| `scripts/plot_validation_results_publication.py` | Generate publication-oriented validation plots. |

## Deprecated or Non-Core Files

`optimize_large_cluster_consensus.py` is retained as an experimental helper and
is not called by `main.nf`. Downstream workflows should not depend on it unless
it is explicitly promoted into the tested pipeline.
