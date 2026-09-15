#!/usr/bin/env nextflow

/*
 * Norovirus wastewater mapping-reference branch.
 *
 * This workflow is independent of the 95%-similarity probe-reference branch.
 * It groups by RdRp_VP1, builds IQ-TREE ML trees with SH-aLRT/UFBoot, performs
 * hybrid ASR/medoid reference selection, and iteratively merges references
 * whose 150 nt windows do not meet the mapping distinguishability threshold.
 */

params.input_file = null
params.output_dir = null
params.threads = 8
params.alignment_min_coverage = 0.5
params.max_n_gap_percent = 10.0
params.iqtree_model = "GTR+F+R4"
params.sh_alrt_replicates = 1000
params.ufboot_replicates = 1000
params.min_sh_alrt = 80.0
params.min_ufboot = 95.0
params.max_within_d95 = 0.05
params.max_within_max = 1.0
params.min_parent_d95 = 0.05
params.min_parent_delta = 0.01
params.asr_pp = 0.90
params.majority_frequency = 0.60
params.max_low_pp_fraction = 0.05
params.max_asr_n_fraction = 0.005
params.max_asr_n_run = 15
params.window_size = 150
params.window_step = 25
params.min_distinguishable_rate = 0.30
params.min_target_mapping_rate = 0.95
params.merge_min_similarity = 0.95
params.max_iterations = 20
params.bowtie2 = "/Users/LuJ/Yunpan/software/Bioinf_software/bowtie2-2.3.3/bowtie2"

if (!params.input_file) error "ERROR: --input_file is required"
if (!params.output_dir) error "ERROR: --output_dir is required"
if (!file(params.input_file).exists()) error "ERROR: input FASTA not found: ${params.input_file}"

def total_cpus = Math.max(1, params.threads as int)
def process_cpus = Math.max(1, Math.min(4, total_cpus) as int)

Channel.fromPath(params.input_file, type: "file").set { input_ch }

process FILTER_INPUT {
    publishDir "${params.output_dir}/01_filtered", mode: "copy"
    cpus 1
    memory "2 GB"

    input:
    path input_fasta

    output:
    path "filtered_input.fasta"
    path "filter_report.txt"

    script:
    """
    python3 ${workflow.projectDir}/scripts/filter_n_sequences.py \
      --input_file ${input_fasta} \
      --output_file filtered_input.fasta \
      --threshold ${params.max_n_gap_percent} > filter_report.txt 2>&1
    """
}

process GROUP_BY_GENOTYPE {
    publishDir "${params.output_dir}/02_groups", mode: "copy"
    cpus 1
    memory "4 GB"

    input:
    path filtered_fasta

    output:
    path "cluster_info.csv"
    path "grouped_sequences.fasta"
    path "genotype_group_report.tsv"

    script:
    """
    python3 ${workflow.projectDir}/scripts/group_sequences_by_rdrp_vp1.py \
      --input ${filtered_fasta} \
      --cluster_info cluster_info.csv \
      --sequences_output grouped_sequences.fasta \
      --report genotype_group_report.tsv
    """
}

process SPLIT_GROUPS {
    publishDir "${params.output_dir}/03_groups", mode: "copy"
    cpus 1
    memory "4 GB"

    input:
    path cluster_info
    path grouped_fasta

    output:
    path "clusters/*"

    script:
    """
    python3 ${workflow.projectDir}/scripts/split_clusters.py \
      --cluster_info ${cluster_info} \
      --sequences ${grouped_fasta} \
      --outdir clusters
    """
}

process ORIENT_AND_ALIGN {
    tag "${group_key}"
    publishDir "${params.output_dir}/04_alignments", mode: "copy", pattern: "aligned_*.fasta"
    cpus process_cpus
    memory "8 GB"
    time "48h"
    maxForks Math.max(1, total_cpus.intdiv(process_cpus))

    input:
    tuple val(group_key), path(group_fasta)

    output:
    tuple val(group_key), path("oriented_*.fasta"), path("aligned_*.fasta")

    script:
    """
    python3 ${workflow.projectDir}/scripts/normalize_orientation.py \
      --input ${group_fasta} \
      --output oriented_${group_fasta.name}

    seq_count=\$(grep -c '^>' oriented_${group_fasta.name})
    if [ "\$seq_count" -eq 1 ]; then
      cp oriented_${group_fasta.name} pretrim_${group_fasta.name}
    else
      python3 -c "from Bio import SeqIO; from Bio.Seq import Seq; r=list(SeqIO.parse('oriented_${group_fasta.name}','fasta')); [setattr(x,'seq',Seq(str(x.seq).replace('-',''))) for x in r]; SeqIO.write(r,'ungapped_${group_fasta.name}','fasta')"
      mafft --quiet --auto --thread ${task.cpus} ungapped_${group_fasta.name} > pretrim_${group_fasta.name}
    fi

    python3 ${workflow.projectDir}/scripts/trim_alignment_ends.py \
      --input pretrim_${group_fasta.name} \
      --output aligned_${group_fasta.name} \
      --min-coverage ${params.alignment_min_coverage}
    """
}

process IQTREE_ASR {
    tag "${group_key}"
    publishDir "${params.output_dir}/05_iqtree", mode: "copy"
    cpus process_cpus
    memory "12 GB"
    time "96h"
    maxForks Math.max(1, total_cpus.intdiv(process_cpus))

    input:
    tuple val(group_key), path(oriented_fasta), path(aligned_fasta)

    output:
    tuple val(group_key), path(oriented_fasta), path(aligned_fasta), path("iqtree_${group_key}")

    script:
    """
    mkdir iqtree_${group_key}
    seq_count=\$(grep -c '^>' ${aligned_fasta})
    if [ "\$seq_count" -ge 4 ]; then
      iqtree3 \
        -s ${aligned_fasta} \
        -m '${params.iqtree_model}' \
        --alrt ${params.sh_alrt_replicates} \
        -B ${params.ufboot_replicates} \
        --bnni \
        --ancestral \
        --asr-min 0.8 \
        -T AUTO \
        --threads-max ${task.cpus} \
        --prefix iqtree_${group_key}/iqtree \
        --redo
    else
      echo "IQ-TREE skipped: \$seq_count sequences" > iqtree_${group_key}/SKIPPED.txt
    fi
    """
}

process SELECT_MAPPING_REFERENCES {
    tag "${group_key}"
    publishDir "${params.output_dir}/06_group_references", mode: "copy"
    cpus process_cpus
    memory "12 GB"
    time "96h"
    maxForks Math.max(1, total_cpus.intdiv(process_cpus))

    input:
    tuple val(group_key), path(oriented_fasta), path(aligned_fasta), path(iqtree_dir)

    output:
    tuple val(group_key), path("mapping_ref_${group_key}")

    script:
    """
    tree_args=""
    if [ -s ${iqtree_dir}/iqtree.treefile ]; then
      tree_args="--tree ${iqtree_dir}/iqtree.treefile --state ${iqtree_dir}/iqtree.state --mldist ${iqtree_dir}/iqtree.mldist"
    fi

    python3 ${workflow.projectDir}/scripts/build_mapping_reference_group.py \
      --alignment ${aligned_fasta} \
      \$tree_args \
      --outdir mapping_ref_${group_key} \
      --bowtie2 ${params.bowtie2} \
      --threads ${task.cpus} \
      --min-sh-alrt ${params.min_sh_alrt} \
      --min-ufboot ${params.min_ufboot} \
      --max-within-d95 ${params.max_within_d95} \
      --max-within-max ${params.max_within_max} \
      --min-parent-d95 ${params.min_parent_d95} \
      --min-parent-delta ${params.min_parent_delta} \
      --asr-pp ${params.asr_pp} \
      --majority-frequency ${params.majority_frequency} \
      --max-low-pp-fraction ${params.max_low_pp_fraction} \
      --max-n-fraction ${params.max_asr_n_fraction} \
      --max-n-run ${params.max_asr_n_run} \
      --window ${params.window_size} \
      --step ${params.window_step} \
      --min-distinguishable-rate ${params.min_distinguishable_rate} \
      --min-target-mapping-rate ${params.min_target_mapping_rate} \
      --merge-min-similarity ${params.merge_min_similarity} \
      --max-iterations ${params.max_iterations}
    """
}

process COLLECT_MAPPING_REFERENCES {
    publishDir "${params.output_dir}/07_final", mode: "copy"
    cpus 1
    memory "4 GB"

    input:
    path group_dirs

    output:
    path "all_mapping_references.fasta"
    path "final_sequence_qc.tsv"
    path "summary.txt"

    script:
    """
    python3 ${workflow.projectDir}/scripts/collect_mapping_references.py \
      --input-glob 'mapping_ref_*/final_mapping_references.fasta' \
      --output all_mapping_references.fasta \
      --report final_sequence_qc.tsv \
      --max-n-gap-percent ${params.max_n_gap_percent}

    {
      echo "Norovirus wastewater mapping-reference branch"
      echo "Input: ${params.input_file}"
      echo "SH-aLRT threshold: ${params.min_sh_alrt}"
      echo "UFBoot threshold: ${params.min_ufboot}"
      echo "Window size: ${params.window_size}"
      echo "Minimum distinguishable-window rate: ${params.min_distinguishable_rate}"
      echo "Final references: \$(grep -c '^>' all_mapping_references.fasta)"
    } > summary.txt
    """
}

workflow {
    FILTER_INPUT(input_ch)
    GROUP_BY_GENOTYPE(FILTER_INPUT.out[0])
    SPLIT_GROUPS(GROUP_BY_GENOTYPE.out[0], GROUP_BY_GENOTYPE.out[1])

    group_ch = SPLIT_GROUPS.out.flatten().map { group_file ->
        tuple(group_file.baseName, group_file)
    }
    aligned_ch = ORIENT_AND_ALIGN(group_ch)
    tree_ch = IQTREE_ASR(aligned_ch)
    reference_ch = SELECT_MAPPING_REFERENCES(tree_ch)
    collected_dirs = reference_ch.map { group_key, group_dir -> group_dir }.collect()
    COLLECT_MAPPING_REFERENCES(collected_dirs)
}
