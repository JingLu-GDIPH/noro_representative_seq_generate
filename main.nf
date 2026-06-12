#!/usr/bin/env nextflow

/*
 * Enhanced Norovirus GII Consensus Sequence Pipeline
 * 
 * IMPORTANT: This pipeline is configured to run all processes in the noro-consensus environment.
 * The conda environment is automatically activated via nextflow.config.
 *
 * This workflow extends the main.nf pipeline by adding iterative consensus validation.
 * After Process 5 (consensus generation), it validates that each consensus file has
 * internal pairwise similarity below 95%. If not, it repeats the consensus process
 * until the threshold is met.
 * 
 * Environment Setup:
 * 1. Create environment: conda env create -f environment.yml
 * 2. Activate environment: conda activate noro-consensus
 * 3. Run pipeline: nextflow run main_consensus_new.nf --input_file <file> --output_dir <dir>
 * 
 * Note: The conda environment is automatically activated for all processes via the config.
 */

// ==================== PARAMETERS ====================
// Required parameters (no default values)
params.input_file = null
params.output_dir = null

// Optional parameters with default values
params.ani_threshold = 0.9
params.cluster_algorithm = "complete"
params.similarity_threshold = 95.0
params.max_iterations = 10  // 最大迭代次数
params.internal_similarity_threshold = 95.0  // 内部相似性阈值
params.alignment_end_min_coverage = 0.5  // 比对两端A/C/G/T最低覆盖率
params.threads = 10

// ==================== PARAMETER VALIDATION ====================
if (!params.input_file) {
    error "ERROR: --input_file parameter is required. Please specify the input FASTA file."
}

if (!params.output_dir) {
    error "ERROR: --output_dir parameter is required. Please specify the output directory."
}

// Validate input file exists
if (!file(params.input_file).exists()) {
    error "ERROR: Input file '${params.input_file}' does not exist."
}

// ==================== LOG INFO ====================
log.info """\
         ENHANCED NOROVIRUS GII CONSENSUS PIPELINE
         ==========================================
         Input File                    : ${params.input_file}
         Output Directory              : ${params.output_dir}
         ANI Threshold                 : ${params.ani_threshold}
         Cluster Algorithm             : ${params.cluster_algorithm}
         Threads                       : ${params.threads}
         Similarity Threshold          : ${params.similarity_threshold}
         Internal Similarity Threshold : ${params.internal_similarity_threshold}
         Alignment End Min Coverage    : ${params.alignment_end_min_coverage}
         Max Iterations                : ${params.max_iterations}
         Enhanced Features             : Iterative consensus validation
         """

// ==================== Input Channel ====================
Channel.fromPath(params.input_file, type: 'file').set { input_fasta_ch }

// ==================== PROCESS 1: Copy and Filter Input FASTA ====================
process MERGE_FASTA {
    tag "Copy and Filter Input FASTA"
    publishDir "${params.output_dir}/01_merged", mode: 'copy'

    input:
    path input_fasta

    output:
    path "merged_sequences.fasta"
    path "filter_report.txt"

    script:
    """
    # 首先复制输入文件
    cp ${input_fasta} input_copy.fasta
    echo "Input FASTA file copied: ${input_fasta}"
    
    # 过滤掉N碱基和gap比例超过10%的序列
    echo "Filtering sequences with N content and gaps > 10%..."
    python3 ${workflow.projectDir}/scripts/filter_n_sequences.py \\
        --input_file input_copy.fasta \\
        --output_file merged_sequences.fasta \\
        --threshold 10.0 > filter_report.txt 2>&1
    
    echo "Sequence filtering completed. Check filter_report.txt for details."
    """
}

// ==================== PROCESS 2: Clustering with vcluster ====================
process CLUSTER {
    tag "vcluster Clustering"
    publishDir "${params.output_dir}/02_cluster", mode: 'copy'

    input:
    path merged_fasta

    output:
    path "cluster_info.csv"
    path "clustered_sequences.fasta"
    path "vcluster_results/*"

    script:
    """
    mkdir -p vcluster_results
    
    echo "步骤1: vcluster预过滤 (使用${task.cpus}线程)..."
    vclust prefilter -i ${merged_fasta} -o vcluster_results/fltr.txt --threads ${task.cpus}
    
    echo "步骤2: 计算成对ANI (使用${task.cpus}线程)..."
    vclust align -i ${merged_fasta} -o vcluster_results/ani.tsv --filter vcluster_results/fltr.txt --threads ${task.cpus}
    
    echo "步骤3: 基于ANI进行聚类..."
    vclust cluster -i vcluster_results/ani.tsv -o vcluster_results/clusters.tsv \
        --ids vcluster_results/ani.ids.tsv --metric ani --ani ${params.ani_threshold} \
        --algorithm ${params.cluster_algorithm}
    
    python3 ${workflow.projectDir}/scripts/process_vcluster_results.py \
        --clusters_file vcluster_results/clusters.tsv \
        --ids_file vcluster_results/ani.ids.tsv \
        --aligned_fasta ${merged_fasta} \
        --cluster_output cluster_info.csv \
        --sequences_output clustered_sequences.fasta
    """
}

// ==================== PROCESS 3: Split Clusters ====================
process SPLIT_CLUSTERS {
    tag "Split Clusters"
    publishDir "${params.output_dir}/03_split", mode: 'copy'

    input:
    path cluster_info
    path clustered_sequences

    output:
    path "clusters/*"

    script:
    """
    python3 ${workflow.projectDir}/scripts/split_clusters.py \
        --cluster_info ${cluster_info} \
        --sequences ${clustered_sequences} \
        --outdir clusters
    """
}

// ==================== PROCESS 4: Multiple Sequence Alignment ====================
process ALIGN {
    tag "MAFFT Alignment"
    publishDir "${params.output_dir}/04_aligned", mode: 'copy', pattern: "aligned_*.fasta"

    input:
    path cluster_fasta

    output:
    tuple path(cluster_fasta), path("aligned_*.fasta")

    script:
    """
    # 获取输入文件名（不包含路径）
    output_name=\$(basename ${cluster_fasta})
    
    seq_count=\$(grep -c "^>" ${cluster_fasta})

    if [ \$seq_count -eq 1 ]; then
        # 单序列不需要多序列比对
        cp ${cluster_fasta} pretrim_\${output_name}
        echo "Single-sequence cluster; alignment skipped for \${output_name}."
    else
        # 等长只表示序列经过padding，不代表同源位点已经对齐。
        # 先移除已有gap，再对每个多序列簇执行MAFFT。
        python3 -c "from Bio import SeqIO; from Bio.Seq import Seq; records=list(SeqIO.parse('${cluster_fasta}', 'fasta')); [setattr(record, 'seq', Seq(str(record.seq).replace('-', ''))) for record in records]; SeqIO.write(records, 'ungapped_\${output_name}', 'fasta')"

        echo "Running MAFFT alignment for \${output_name} (\$seq_count sequences)..."
        mafft --thread ${task.cpus} --auto ungapped_\${output_name} > pretrim_\${output_name}
    fi

    # 仅从比对两端删除A/C/G/T覆盖率低于阈值的连续列。
    # N和gap不计为有效覆盖，内部低覆盖列保留。
    python3 ${workflow.projectDir}/scripts/trim_alignment_ends.py \\
        --input pretrim_\${output_name} \\
        --output aligned_\${output_name} \\
        --min-coverage ${params.alignment_end_min_coverage}

    # 验证修剪后的比对结果
    python3 -c "from Bio import SeqIO; records=list(SeqIO.parse('aligned_\${output_name}', 'fasta')); lengths={len(record.seq) for record in records}; assert len(records)==\$seq_count, 'Alignment sequence count changed'; assert len(lengths)==1, 'Trimmed output is not aligned'; assert next(iter(lengths)) > 0, 'Trimmed alignment is empty'; print(f'Alignment completed: {len(records)} sequences, trimmed alignment length {next(iter(lengths))}')"
    """
}

// ==================== PROCESS 5: Generate Consensus Sequences (Molecular Phylogeny) ====================
process CONSENSUS {
    tag "Consensus Generation"
    publishDir "${params.output_dir}/05_consensus", mode: 'copy'

    input:
    tuple path(cluster_fasta), path(aligned_fasta)

    output:
    path "*_consensus_*.fasta"

    script:
    """
    echo "Consensus logic revision: 2026-06-12-mafft-state-root-v2"

    # 获取文件名信息
    filename=\$(basename ${cluster_fasta})
    basename=\${filename%.*}
    
    # 检查序列数量
    sequence_count=\$(grep -c "^>" ${cluster_fasta})
    
    # 从文件名提取基因型和cluster信息
    # 处理格式: cluster_cluster_0_GI.6.fasta
    if [[ \$filename =~ cluster_cluster_([0-9]+)_([^.]+)\\.fasta ]]; then
        cluster_id=\${BASH_REMATCH[1]}
        genotype=\${BASH_REMATCH[2]}
    else
        # 使用sed进行更精确的解析
        cluster_id=\$(echo \$filename | sed -n 's/.*cluster_cluster_\\([0-9]*\\)_.*\\.fasta/\\1/p')
        genotype=\$(echo \$filename | sed -n 's/.*cluster_cluster_[0-9]*_\\([A-Z0-9.]*\\)\\.fasta/\\1/p')
        
        # 如果解析失败，使用默认值
        if [ -z "\$cluster_id" ] || [ "\$cluster_id" = "\$filename" ]; then
            cluster_id="unknown"
        fi
        if [ -z "\$genotype" ] || [ "\$genotype" = "\$filename" ]; then
            genotype="unknown"
        fi
    fi
    
    echo "解析结果: cluster_id=\$cluster_id, genotype=\$genotype"
    
    # 检查基因型是否包含"mixed"，如果是则保持"mixed"命名
    # 注意：不要将包含下划线的基因型都重命名为"mixed"
    # 因为enhanced_consensus_with_validation.py会生成包含下划线的正确基因型名称
    if [[ \$genotype == *mixed* ]]; then
        genotype="mixed"
    fi
    
    echo "最终基因型: \$genotype"
    
    consensus_file="\${genotype}_consensus_\${cluster_id}.fasta"
    
    if [ "\$sequence_count" -eq 1 ]; then
        cp ${aligned_fasta} \$consensus_file
    else
        # 使用generate_consensus.py脚本进行分子进化树分析：
        # 1. 使用IQ-TREE构建系统发育树和祖先序列重建
        # 2. 进行midpoint rooting和节点编号分配
        # 3. 寻找optimal node（满足相似性阈值的节点）
        # 4. 基于optimal node生成一致性序列
        python3 ${workflow.projectDir}/scripts/generate_consensus.py \\
            --aligned_fasta ${aligned_fasta} \\
            --consensus_output \$consensus_file \\
            --similarity_threshold ${params.similarity_threshold} \\
            --threads ${task.cpus}
    fi
    """
}

// ==================== PROCESS 6: Enhanced Consensus Generation with Validation (Molecular Phylogeny) ====================
process ENHANCED_CONSENSUS {
    tag "Enhanced Consensus Generation with Validation"
    publishDir "${params.output_dir}/06_enhanced_consensus", mode: 'copy'

    input:
    path consensus_file

    output:
    path "final_consensus_*.fasta"

    script:
    """
    echo "Enhanced consensus logic revision: 2026-06-12-redundancy-pruning-v5"

    # 获取文件名信息
    filename=\$(basename ${consensus_file})
    basename=\${filename%.*}
    
    echo "Enhanced consensus generation with validation for \$filename..."
    
    # 检查序列数量
    sequence_count=\$(grep -c "^>" ${consensus_file})
    
    if [ "\$sequence_count" -le 1 ]; then
        # Single sequence, no processing needed, just copy the file
        cp ${consensus_file} final_consensus_\${basename}.fasta
    else
        # 使用修正后的enhanced_consensus_with_validation.py脚本
        # 该脚本实现：
        # 1. 迭代验证Process 5输出的内部相似性
        # 2. 如果相似性>95%，调用分子进化树构建流程
        # 3. 使用IQ-TREE构建系统发育树和祖先序列重建
        # 4. 寻找optimal node并生成一致性序列
        # 5. 重复直到所有序列相似性<=95%
        python3 ${workflow.projectDir}/scripts/enhanced_consensus_with_validation.py \\
            --input_file ${consensus_file} \\
            --output_file final_consensus_\${basename}.fasta \\
            --similarity_threshold ${params.internal_similarity_threshold} \\
            --max_iterations ${params.max_iterations} \\
            --threads ${task.cpus}
    fi

    # 最终代表序列不保留由比对产生的两端N或gap；内部位点保持不变。
    python3 ${workflow.projectDir}/scripts/trim_sequence_ends.py \\
        --input final_consensus_\${basename}.fasta \\
        --output final_consensus_\${basename}.fasta
    """
}

// ==================== PROCESS 7: Collect All Results ====================
process COLLECT_RESULTS {
    tag "Collect All Results"
    publishDir "${params.output_dir}/07_final_results", mode: 'copy'

    input:
    path fasta_files

    output:
    path "all_final_consensus.fasta"
    path "summary_report.txt"

    script:
    """
    #!/bin/bash
    set -e
    
    echo "Collecting all final consensus sequences..."
    
    # Combine all consensus files and keep FASTA identifiers globally unique.
    python3 ${workflow.projectDir}/scripts/combine_consensus_fastas.py \\
        --input_glob "final_consensus_*.fasta" \\
        --output_file all_final_consensus.raw.fasta

    # Remove near-identical representatives that remain across independent clusters.
    python3 ${workflow.projectDir}/scripts/prune_redundant_consensus.py \\
        --input all_final_consensus.raw.fasta \\
        --output all_final_consensus.fasta \\
        --threshold ${params.internal_similarity_threshold} \\
        --threads ${task.cpus}
    
    # Generate summary
    echo "Enhanced Norovirus GII Consensus Pipeline - Final Summary" > summary_report.txt
    echo "=========================================================" >> summary_report.txt
    echo "Date: \$(date)" >> summary_report.txt
    echo "Input file: ${params.input_file}" >> summary_report.txt
    echo "Internal similarity threshold: ${params.internal_similarity_threshold}%" >> summary_report.txt
    echo "Max iterations: ${params.max_iterations}" >> summary_report.txt
    echo "" >> summary_report.txt
    echo "Final consensus files: \$(ls final_consensus_*.fasta | wc -l)" >> summary_report.txt
    echo "Total final sequences: \$(grep -c '^>' all_final_consensus.fasta)" >> summary_report.txt
    """
}

// ==================== PROCESS 8: Final Validation ====================
process FINAL_VALIDATION {
    tag "Final Validation of All Consensus Sequences"
    publishDir "${params.output_dir}/08_validation", mode: 'copy'

    input:
    path "all_final_consensus.fasta"

    output:
    path "final_validation_report.txt"

    script:
    """
    #!/bin/bash
    set -e
    
    echo "Performing final validation of all consensus sequences..."
    echo "Equal-length alignments will be validated directly."
    
    # Run final validation
    python3 ${workflow.projectDir}/scripts/validate_consensus_internal.py \\
        --consensus_file all_final_consensus.fasta \\
        --output_report final_validation_report.txt \\
        --similarity_threshold ${params.internal_similarity_threshold}
    
    echo "Final validation completed"
    """
}

// ==================== Workflow Definition ====================
workflow {
    // 原始流程
    MERGE_FASTA(input_fasta_ch)
    CLUSTER(MERGE_FASTA.out[0])
    SPLIT_CLUSTERS(CLUSTER.out[0], CLUSTER.out[1])
    
    // 初始比对和一致性序列生成
    initial_align_ch = SPLIT_CLUSTERS.out.flatten() | ALIGN
    initial_consensus_ch = initial_align_ch | CONSENSUS
    
    // 对process 5的输出进行增强一致性序列生成和验证
    final_consensus_ch = initial_consensus_ch | ENHANCED_CONSENSUS
    
    // 收集所有最终结果 - 参照continue_pipeline.nf的语法
    all_fasta = final_consensus_ch.map { it }.collect()
    all_consensus_ch = COLLECT_RESULTS(all_fasta)
    
    // 最终验证
    validation_ch = FINAL_VALIDATION(all_consensus_ch[0])
} 
