#!/usr/bin/env python3
"""
Enhanced Consensus Generation with Validation
This script implements Process 6 logic:
1. First validate pairwise similarity of Process 5 output
2. If minimum similarity > 95%, perform new round of consensus generation
3. Loop until all sequences have similarity <= 95%
4. Use JC model and 80% posterior probability threshold for IQ-TREE
"""

import argparse
import csv
import subprocess
import os
import sys
import tempfile
import logging
import re
import shutil
from collections import defaultdict, Counter
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import multiprocessing as mp
from functools import partial

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

VALID_BASES = set("ACGT")


def make_record_ids_unique(records):
    """Return records with stable, unique FASTA identifiers."""
    seen = Counter()
    unique_records = []

    for record in records:
        original_id = record.id
        seen[original_id] += 1
        if seen[original_id] > 1:
            record.id = f"{original_id}_dup{seen[original_id]}"
            record.name = record.id
            record.description = ""
        unique_records.append(record)

    return unique_records


def normalize_fasta_ids(fasta_file):
    """Rewrite a FASTA file in place so IQ-TREE always receives unique IDs."""
    records = make_record_ids_unique(list(SeqIO.parse(fasta_file, "fasta")))
    temp_file = f"{fasta_file}.unique.tmp"
    SeqIO.write(records, temp_file, "fasta")
    os.replace(temp_file, fasta_file)


def calculate_sequence_similarity(seq1, seq2):
    """
    计算两条序列的相似性百分比
    只有当两个位置都为 A/C/G/T 时才计入分母；N 和 gap 均不作为有效覆盖。
    
    Args:
        seq1 (str): 第一条序列
        seq2 (str): 第二条序列
        
    Returns:
        float: 相似性百分比
    """
    if len(seq1) != len(seq2):
        raise ValueError("序列长度不匹配")
    
    matches = 0
    total_positions = 0
    
    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        if base1 not in VALID_BASES or base2 not in VALID_BASES:
            continue
        total_positions += 1
        if base1 == base2:
            matches += 1
    
    if total_positions == 0:
        return 0.0
    
    similarity = (matches / total_positions) * 100
    return similarity

def calculate_similarity_wrapper(args):
    """
    包装函数，用于多进程计算相似性
    """
    i, j, seq1, seq2 = args
    return calculate_sequence_similarity(seq1, seq2)

def calculate_pairwise_similarity_parallel(sequences, num_processes=None):
    """
    使用多进程并行计算所有序列对的相似性
    
    Args:
        sequences (list): 序列列表
        num_processes (int): 进程数，默认为CPU核心数
        
    Returns:
        list: 所有相似性值
    """
    if num_processes is None:
        num_processes = min(mp.cpu_count(), 8)  # 限制最大进程数
    
    if len(sequences) <= 1:
        return [100.0]
    
    # 生成所有序列对
    pairs = []
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            pairs.append((i, j, sequences[i], sequences[j]))
    
    if not pairs:
        return [100.0]
    
    # 使用进程池并行计算
    with mp.Pool(processes=num_processes) as pool:
        similarities = pool.map(calculate_similarity_wrapper, pairs)
    
    return similarities

def calculate_maximum_pairwise_similarity(fasta_file):
    """
    计算FASTA文件中所有序列的最高pairwise相似性（优化版本）
    
    Args:
        fasta_file (str): 输入FASTA文件路径
        
    Returns:
        float: 最高相似性百分比
    """
    try:
        # 读取比对文件
        alignment = AlignIO.read(fasta_file, "fasta")
        
        if len(alignment) <= 1:
            return 100.0  # 单个序列，相似性为100%
        
        logging.info(f"开始计算 {len(alignment)} 条序列的最高pairwise相似性...")
        
        # 提取序列字符串
        sequences = [str(record.seq) for record in alignment]
        
        # 使用并行计算
        similarities = calculate_pairwise_similarity_parallel(sequences)
        
        if not similarities:
            return 0.0
        
        max_similarity = max(similarities)
        logging.info(f"计算得到最高pairwise相似性: {max_similarity:.2f}%")
        
        return max_similarity
        
    except Exception as e:
        logging.error(f"计算pairwise相似性时发生错误: {str(e)}")
        import traceback
        logging.error(f"错误详情: {traceback.format_exc()}")
        return 0.0

def run_iqtree_enhanced(aligned_file, output_prefix, threads=8, fast_mode=False):
    """
    运行IQ-TREE进行系统发育树重建和祖先序列重建
    """
    try:
        # 根据模式选择bootstrap次数
        bootstrap_count = 50 if fast_mode else 100
        
        # 运行IQ-TREE（优化版本，减少bootstrap次数）
        cmd = [
            'iqtree', '-s', aligned_file,
            '-m', 'JC',  # 使用JC核酸替换模型
            '--ancestral',  # 生成祖先序列
            '--asr-min', '0.8',  # 后验概率阈值为80%
            '-nt', str(threads),
            '-pre', output_prefix,
            '-quiet'
        ]
        
        logging.info(f"运行IQ-TREE命令: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # 检查输出文件
        tree_file = f"{output_prefix}.treefile"
        ancestral_file = f"{output_prefix}.state"
        
        if not os.path.exists(tree_file):
            raise FileNotFoundError(f"IQ-TREE树文件未生成: {tree_file}")
        
        if not os.path.exists(ancestral_file):
            raise FileNotFoundError(f"IQ-TREE祖先序列文件未生成: {ancestral_file}")
        
        logging.info(f"IQ-TREE运行成功，生成文件: {tree_file}, {ancestral_file}")
        return tree_file, ancestral_file
        
    except subprocess.CalledProcessError as e:
        logging.error(f"IQ-TREE运行失败: {e}")
        logging.error(f"错误输出: {e.stderr}")
        raise
    except Exception as e:
        logging.error(f"运行IQ-TREE时发生错误: {str(e)}")
        raise

def midpoint_root_tree_manual(tree_file, output_file):
    """
    手动实现midpoint rooting
    """
    try:
        # 读取树文件
        tree = Phylo.read(tree_file, 'newick')
        
        # 找到最长的路径
        terminals = list(tree.get_terminals())
        max_distance = 0
        longest_path = None
        
        for i, term1 in enumerate(terminals):
            for j, term2 in enumerate(terminals[i+1:], i+1):
                distance = tree.distance(term1, term2)
                if distance > max_distance:
                    max_distance = distance
                    longest_path = (term1, term2)
        
        if longest_path:
            # 在最长路径的中点进行rooting
            midpoint_distance = max_distance / 2
            logging.info(f"最长路径: {longest_path[0].name} -> {longest_path[1].name}, 距离: {max_distance:.6f}")
            logging.info(f"中点距离: {midpoint_distance:.6f}")
            
            # 重新root树
            tree.root_at_midpoint()
            
            # 保存rooted树
            Phylo.write(tree, output_file, 'newick')
            logging.info(f"Midpoint rooting完成，结果保存在: {output_file}")
        else:
            # 如果无法找到最长路径，直接复制原文件
            import shutil
            shutil.copy2(tree_file, output_file)
            logging.warning(f"无法找到最长路径，使用原始树文件: {output_file}")
        
        return output_file
        
    except Exception as e:
        logging.error(f"Midpoint rooting失败: {str(e)}")
        # 如果失败，直接复制原文件
        import shutil
        shutil.copy2(tree_file, output_file)
        logging.warning(f"使用原始树文件: {output_file}")
        return output_file

def assign_node_numbers(tree):
    """
    为树中的节点分配编号
    """
    node_counter = 1
    
    # 为内部节点分配编号
    for node in tree.get_nonterminals():
        if not hasattr(node, 'name') or not node.name:
            node.name = f"Node_{node_counter}"
            node_counter += 1
    
    # 为末端节点分配编号（如果没有名称）
    for node in tree.get_terminals():
        if not hasattr(node, 'name') or not node.name:
            node.name = f"Tip_{node_counter}"
            node_counter += 1
    
    logging.info(f"为 {node_counter-1} 个节点分配了编号")
    return tree

def parse_ancestral_sequences(ancestral_file):
    """Parse IQ-TREE's tab-separated .state ancestral reconstruction file."""
    states_by_node = defaultdict(dict)

    try:
        with open(ancestral_file, "r") as handle:
            data_lines = (line for line in handle if not line.startswith("#"))
            reader = csv.DictReader(data_lines, delimiter="\t")
            for row in reader:
                states_by_node[row["Node"]][int(row["Site"])] = row["State"]

        ancestral_sequences = {
            node: "".join(states[site] for site in sorted(states))
            for node, states in states_by_node.items()
        }
        logging.info(f"解析了 {len(ancestral_sequences)} 个祖先序列")
        return ancestral_sequences
    except Exception as e:
        logging.error(f"解析祖先序列文件失败: {str(e)}")
        return {}

def parse_tip_sequences(aligned_file):
    """
    解析比对文件中的末端序列
    """
    tip_sequences = {}
    
    try:
        alignment = AlignIO.read(aligned_file, "fasta")
        
        for record in alignment:
            tip_sequences[record.id] = str(record.seq)
        
        logging.info(f"解析了 {len(tip_sequences)} 个末端序列")
        return tip_sequences
        
    except Exception as e:
        logging.error(f"解析末端序列失败: {str(e)}")
        return {}


def extract_norovirus_genotype(text):
    """Extract a combined RdRp_VP1 label, or fall back to one genotype label."""
    pair_match = re.search(
        r"(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)",
        str(text),
    )
    if pair_match:
        return f"{pair_match.group(1)}_{pair_match.group(2)}"

    match = re.search(r'(?<![A-Z0-9])G(?:I|II|IX)\.P?[A-Za-z0-9]+(?![A-Za-z0-9])', str(text))
    return match.group() if match else None


def extract_common_genotype_from_cluster(aligned_file):
    """
    从比对文件名中提取基因型信息
    """
    try:
        # 从文件名中提取基因型
        filename = os.path.basename(aligned_file)
        
        genotype = extract_norovirus_genotype(filename)
        if genotype:
            return genotype
        
        # 如果文件名中没有基因型信息，从序列ID中提取
        alignment = AlignIO.read(aligned_file, "fasta")
        genotypes = []
        
        for record in alignment:
            genotype = extract_norovirus_genotype(record.id)
            if genotype:
                genotypes.append(genotype)
        
        if genotypes:
            # 返回最常见的基因型
            most_common = Counter(genotypes).most_common(1)[0][0]
            return most_common
        
        return "unknown"
        
    except Exception as e:
        logging.error(f"提取基因型信息失败: {str(e)}")
        return "unknown"

def extract_cluster_number_from_filename(filename):
    """从文件名中提取cluster序号，兼容多种格式"""
    import re
    basename = os.path.basename(filename)
    # 1. GI.3_consensus_1.fasta 或 GI.3_consensus_1
    m = re.search(r'_consensus_([^.]+)', basename)
    if m:
        return m.group(1)
    # 2. cluster_cluster_1_GI.3.fasta
    m = re.search(r'cluster_cluster_([0-9]+)_', basename)
    if m:
        return m.group(1)
    # 3. GI.3_cluster_1_node_1.fasta
    m = re.search(r'_cluster_([0-9]+)', basename)
    if m:
        return m.group(1)
    return "unknown"

def extract_enhanced_iteration_from_filename(filename):
    """
    从文件名中提取enhanced迭代次数
    """
    try:
        # 尝试从文件名中提取enhanced迭代次数
        match = re.search(r'enhanced(\d+)', filename)
        if match:
            return int(match.group(1))
        
        # 如果没有找到，返回1（第一次迭代）
        return 1
        
    except Exception as e:
        logging.error(f"提取enhanced迭代次数失败: {str(e)}")
        return 1

def get_node_genotypes(tree, node, ancestral_sequences, tip_sequences):
    """
    获取节点包含的所有基因型信息
    """
    genotypes = set()
    
    try:
        # 如果是末端节点
        if node.is_terminal():
            if node.name in tip_sequences:
                genotype = extract_norovirus_genotype(node.name)
                if genotype:
                    genotypes.add(genotype)
        else:
            # 对于内部节点，获取所有后代节点的基因型
            for tip in tree.get_terminals():
                # 检查tip是否是node的后代
                is_descendant = False
                current = tip
                while current != tree.root:
                    # 找到current的父节点
                    parent = None
                    for potential_parent in tree.get_nonterminals():
                        if potential_parent != current:
                            for child in potential_parent.clades:
                                if child == current:
                                    parent = potential_parent
                                    break
                            if parent:
                                break
                    
                    if parent == node:
                        is_descendant = True
                        break
                    elif parent is None:
                        break
                    current = parent
                
                if is_descendant and tip.name in tip_sequences:
                    genotype = extract_norovirus_genotype(tip.name)
                    if genotype:
                        genotypes.add(genotype)
        
        return list(genotypes)
        
    except Exception as e:
        logging.error(f"获取节点基因型失败: {str(e)}")
        return []

def get_node_genotype_name(genotypes):
    """
    根据基因型列表生成节点基因型名称
    """
    if not genotypes:
        return "unknown"
    
    if len(genotypes) == 1:
        # 如果只有一个基因型，直接返回该基因型名称
        return genotypes[0]
    
    # 如果有多个基因型，使用mixed格式
    return "mixed_" + str(len(genotypes)) + "_genotypes"

def generate_simple_consensus_from_sequences(sequences):
    """
    从序列列表生成简单的一致性序列
    """
    if not sequences:
        return ""
    
    if len(sequences) == 1:
        return sequences[0]
    
    # 如果是两条序列，检查相似性
    if len(sequences) == 2:
        similarity = calculate_sequence_similarity(sequences[0], sequences[1])
        logging.info(f"两条序列的相似性: {similarity:.2f}%")
        
        # 如果相似性>=95%，选择第一条序列作为代表
        if similarity >= 95.0:
            logging.info(f"两条序列相似性 >= 95%，选择第一条序列作为代表")
            return sequences[0]
        else:
            logging.info(f"两条序列相似性 < 95%，生成一致性序列")
    
    # 获取序列长度
    seq_length = len(sequences[0])
    
    # 逐位计算一致性
    consensus = ""
    for pos in range(seq_length):
        # 收集该位置的所有碱基
        bases = []
        for seq in sequences:
            if pos < len(seq):
                base = seq[pos].upper()
                if base in VALID_BASES:
                    bases.append(base)
        
        if not bases:
            consensus += '-'
        else:
            # 计算最常见的碱基
            base_counts = Counter(bases)
            most_common_base = base_counts.most_common(1)[0][0]
            consensus += most_common_base
    
    return consensus


def generate_majority_consensus_from_aligned_sequences(sequences):
    """Generate a majority consensus without choosing one record as the winner."""
    if not sequences:
        return ""

    seq_length = max(len(sequence) for sequence in sequences)
    consensus = []
    for pos in range(seq_length):
        bases = []
        for sequence in sequences:
            if pos < len(sequence):
                base = sequence[pos].upper()
                if base in VALID_BASES:
                    bases.append(base)
        if not bases:
            consensus.append("-")
        else:
            consensus.append(Counter(bases).most_common(1)[0][0])
    return "".join(consensus)


def align_records_for_similarity(records, threads=1):
    """Return aligned sequence strings in input order for similarity comparisons."""
    if len(records) <= 1:
        return [str(record.seq).upper() for record in records], "not needed"

    lengths = {len(record.seq) for record in records}
    if len(lengths) == 1:
        return [str(record.seq).upper() for record in records], "existing"

    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = os.path.join(temp_dir, "input.fasta")
        output_path = os.path.join(temp_dir, "aligned.fasta")
        SeqIO.write(records, input_path, "fasta")

        with open(output_path, "w") as output_handle:
            result = subprocess.run(
                ["mafft", "--quiet", "--auto", "--thread", str(threads), input_path],
                stdout=output_handle,
                stderr=subprocess.PIPE,
                text=True,
            )

        if result.returncode != 0:
            raise RuntimeError(
                f"MAFFT failed while aligning consensus records: {result.stderr[:500]}"
            )

        aligned_records = list(SeqIO.parse(output_path, "fasta"))
        aligned_by_id = {record.id: str(record.seq).upper() for record in aligned_records}
        missing = [record.id for record in records if record.id not in aligned_by_id]
        if missing:
            raise RuntimeError(f"MAFFT dropped sequence IDs: {missing[:5]}")
        return [aligned_by_id[record.id] for record in records], "MAFFT"


def write_aligned_records(records, aligned_sequences, output_path):
    """Write aligned sequences using the original record IDs."""
    aligned_records = []
    for record, sequence in zip(records, aligned_sequences):
        aligned_records.append(SeqRecord(Seq(sequence), id=record.id, description=""))
    SeqIO.write(aligned_records, output_path, "fasta")


def parse_cluster_number_from_id(sequence_id):
    match = re.search(r"cluster[_-]([0-9]+)", sequence_id)
    return match.group(1) if match else None


def merge_consensus_id(records, iteration, merge_index):
    """Build a stable ID for a merged consensus representative."""
    cluster_numbers = [
        parse_cluster_number_from_id(record.id)
        for record in records
        if parse_cluster_number_from_id(record.id)
    ]
    cluster_number = Counter(cluster_numbers).most_common(1)[0][0] if cluster_numbers else "NA"

    first_id = records[0].id
    prefix_match = re.match(r"(.+?)_cluster[_-][0-9]+", first_id)
    if prefix_match:
        prefix = prefix_match.group(1)
    else:
        labels = [extract_norovirus_genotype(record.id) for record in records]
        labels = [label for label in labels if label]
        prefix = Counter(labels).most_common(1)[0][0] if labels else "unknown"

    return f"{prefix}_cluster_{cluster_number}_enhanced{iteration}_merge_{merge_index}"


def find_high_similarity_pairs(aligned_sequences, threshold):
    """Return pairs whose pairwise similarity remains at or above threshold."""
    pairs = []
    for i in range(len(aligned_sequences)):
        for j in range(i + 1, len(aligned_sequences)):
            similarity = calculate_sequence_similarity(aligned_sequences[i], aligned_sequences[j])
            if similarity >= threshold:
                pairs.append((similarity, i, j))
    return pairs


def tree_distances_for_pairs(records, aligned_sequences, pairs, threads):
    """Estimate tree distances for high-similarity pairs; fall back to zero on failure."""
    if not pairs:
        return {}

    with tempfile.TemporaryDirectory() as temp_dir:
        aligned_path = os.path.join(temp_dir, "merge_candidates.aligned.fasta")
        prefix = os.path.join(temp_dir, "merge_tree")
        write_aligned_records(records, aligned_sequences, aligned_path)
        try:
            tree_file, _ = run_iqtree_enhanced(aligned_path, prefix, threads, fast_mode=True)
            tree = Phylo.read(tree_file, "newick")
            distances = {}
            for _, i, j in pairs:
                try:
                    distances[(i, j)] = tree.distance(records[i].id, records[j].id)
                except Exception:
                    distances[(i, j)] = 0.0
            return distances
        except Exception as exc:
            logging.warning(f"树引导相似consensus合并建树失败，回退到相似性排序: {exc}")
            return {(i, j): 0.0 for _, i, j in pairs}


def tree_guided_merge_similar_consensus(fasta_file, threshold=95.0, threads=1, iteration=1):
    """
    Merge highly similar consensus representatives instead of deleting one.

    The previous pruning step removed later records greedily. This routine keeps
    all information by replacing each selected high-similarity pair with a newly
    rebuilt majority consensus. Pairs are prioritized by sequence similarity and,
    when possible, by short patristic distance in a consensus-only tree.
    """
    records = make_record_ids_unique(list(SeqIO.parse(fasta_file, "fasta")))
    if len(records) <= 1:
        return 0

    aligned_sequences, alignment_method = align_records_for_similarity(records, threads)
    pairs = find_high_similarity_pairs(aligned_sequences, threshold)
    if not pairs:
        logging.info("树引导合并检查: 未发现高于阈值的consensus相似对")
        return 0

    tree_distances = tree_distances_for_pairs(records, aligned_sequences, pairs, threads)
    pairs.sort(key=lambda item: (-item[0], tree_distances.get((item[1], item[2]), 0.0)))

    used = set()
    merge_pairs = []
    for similarity, i, j in pairs:
        if i in used or j in used:
            continue
        merge_pairs.append((similarity, i, j))
        used.add(i)
        used.add(j)

    if not merge_pairs:
        return 0

    merged_records = []
    merge_lookup = {}
    for merge_index, (similarity, i, j) in enumerate(merge_pairs, start=1):
        merge_id = merge_consensus_id([records[i], records[j]], iteration, merge_index)
        sequence = generate_majority_consensus_from_aligned_sequences(
            [aligned_sequences[i], aligned_sequences[j]]
        )
        merged_records.append(SeqRecord(Seq(sequence), id=merge_id, description=""))
        merge_lookup[i] = merge_id
        merge_lookup[j] = merge_id
        logging.info(
            "树引导合并相似consensus: %s + %s -> %s (similarity %.2f%%, tree distance %.6f)",
            records[i].id,
            records[j].id,
            merge_id,
            similarity,
            tree_distances.get((i, j), 0.0),
        )

    output_records = []
    emitted_merges = set()
    for index, record in enumerate(records):
        if index in merge_lookup:
            merge_id = merge_lookup[index]
            if merge_id not in emitted_merges:
                output_records.append(next(record for record in merged_records if record.id == merge_id))
                emitted_merges.add(merge_id)
        else:
            output_records.append(record)

    temp_file = f"{fasta_file}.tree_merge.tmp"
    SeqIO.write(make_record_ids_unique(output_records), temp_file, "fasta")
    shutil.move(temp_file, fasta_file)
    logging.info(
        "树引导合并完成: %d 条输入代表序列 -> %d 条输出代表序列；合并 %d 对；比对方法: %s",
        len(records),
        len(output_records),
        len(merge_pairs),
        alignment_method,
    )
    return len(merge_pairs)

def generate_enhanced_consensus_with_node_analysis(aligned_file, consensus_file, iteration, similarity_threshold=95.0, threads=8):
    """
    使用节点分析生成enhanced consensus序列
    参考generate_consensus.py的generate_consensus_with_node_analysis方法
    """
    try:
        logging.info(f"开始生成enhanced consensus，迭代次数: {iteration}")
        
        # 创建临时目录
        with tempfile.TemporaryDirectory() as temp_dir:
            # 步骤1: 运行IQ-TREE
            output_prefix = os.path.join(temp_dir, "iqtree_output")
            tree_file, ancestral_file = run_iqtree_enhanced(aligned_file, output_prefix, threads)
            
            # 步骤2: Midpoint rooting
            rooted_tree_file = os.path.join(temp_dir, "rooted_tree.treefile")
            midpoint_root_tree_manual(tree_file, rooted_tree_file)
            
            # 步骤3: 分配节点编号
            tree = Phylo.read(rooted_tree_file, 'newick')
            tree = assign_node_numbers(tree)
            
            # 步骤4: 解析序列
            ancestral_sequences = parse_ancestral_sequences(ancestral_file)
            tip_sequences = parse_tip_sequences(aligned_file)
            
            # 步骤5: 提取基因型和cluster信息
            genotype = extract_common_genotype_from_cluster(aligned_file)
            cluster_number = extract_cluster_number_from_filename(os.path.basename(aligned_file))
            
            logging.info(f"基因型: {genotype}, Cluster: {cluster_number}")
            
            # 步骤6: 寻找最优节点
            optimal_nodes = find_optimal_nodes(tree, ancestral_sequences, tip_sequences, similarity_threshold)
            
            if not optimal_nodes:
                logging.warning("未找到满足条件的节点，将保留所有原始序列")
                # 直接复制输入文件
                with open(aligned_file, 'r') as src, open(consensus_file, 'w') as dst:
                    dst.write(src.read())
                return
            
            # 步骤7: 生成consensus序列
            consensus_sequences = []
            for i, (node, node_similarity, is_terminal) in enumerate(optimal_nodes):
                logging.info(f"处理节点 {i+1}: {node.name}, 节点最小相似性: {node_similarity:.2f}%, 是否为末端节点: {is_terminal}")
                
                # 获取节点的序列
                if node.name in ancestral_sequences:
                    optimal_sequence = ancestral_sequences[node.name]
                else:
                    # 如果没有祖先序列，使用该节点包含的所有序列生成一致性序列
                    node_seqs = get_node_sequences(tree, node, ancestral_sequences, tip_sequences)
                    if node_seqs:
                        optimal_sequence = generate_simple_consensus_from_sequences(node_seqs)
                    else:
                        logging.warning(f"无法获取节点 {node.name} 的序列，跳过")
                        continue
                
                # 获取节点内部的基因型信息
                node_genotypes = get_node_genotypes(tree, node, ancestral_sequences, tip_sequences)
                node_genotype_name = get_node_genotype_name(node_genotypes)
                logging.info(f"节点 {node.name} 包含的基因型: {node_genotypes}, 节点基因型名称: {node_genotype_name}")
                
                # 生成consensus序列ID
                node_num = i + 1
                
                # 检查是否是末端节点且不满足阈值条件（直接输出末端序列）
                should_use_original_name = False
                if is_terminal:
                    # 检查该末端节点的父节点是否满足阈值
                    parent_satisfies_threshold = False
                    for potential_parent in tree.get_nonterminals():
                        if potential_parent != node:
                            for child in potential_parent.clades:
                                if child == node:
                                    # 找到父节点，检查父节点是否满足阈值
                                    parent_seqs = get_node_sequences(tree, potential_parent, ancestral_sequences, tip_sequences)
                                    if parent_seqs:
                                        parent_similarity = calculate_node_similarity(parent_seqs)
                                        if parent_similarity >= similarity_threshold:
                                            parent_satisfies_threshold = True
                                            break
                            if parent_satisfies_threshold:
                                break
                    
                    # 如果父节点不满足阈值，则使用原始名称
                    if not parent_satisfies_threshold:
                        should_use_original_name = True
                
                if should_use_original_name:
                    # 对于不满足条件的末端节点，保留原来的序列名称
                    if node.name in tip_sequences:
                        consensus_id = node.name
                    else:
                        # 如果无法获取原始名称，使用默认格式
                        consensus_id = f"{node_genotype_name}_cluster_{cluster_number}_enhanced{iteration}_node_{node_num}"
                else:
                    # 对于满足条件的节点或内部节点，使用新的命名格式
                    # 使用从文件名中提取的enhanced迭代次数，而不是传入的iteration参数
                    consensus_id = f"{node_genotype_name}_cluster_{cluster_number}_enhanced{iteration}_node_{node_num}"
                
                consensus_sequences.append((consensus_id, optimal_sequence))
                logging.info(f"生成consensus序列: {consensus_id}")
            
            # 保存所有consensus序列
            with open(consensus_file, "w") as f:
                for consensus_id, sequence in consensus_sequences:
                    f.write(f">{consensus_id}\n{sequence}\n")
            
            logging.info(f"共生成 {len(consensus_sequences)} 个consensus序列，已保存至: {consensus_file}")
    
    except Exception as e:
        logging.error(f"生成enhanced consensus失败: {str(e)}")
        import traceback
        logging.error(f"错误详情: {traceback.format_exc()}")
        raise

def get_node_sequences(tree, node, ancestral_sequences, tip_sequences):
    """
    获取指定节点包含的所有序列（包括祖先序列和后代序列）
    参考generate_consensus.py的实现
    """
    sequences = []
    
    # 如果是末端节点，直接返回该节点的序列
    if node.is_terminal():
        if node.name in tip_sequences:
            sequences.append(tip_sequences[node.name])
        return sequences
    
    # 对于内部节点，获取该节点的所有后代节点
    descendants = []
    for tip in tree.get_terminals():
        # 检查tip是否是node的后代
        is_descendant = False
        current = tip
        while current != tree.root:
            # 找到current的父节点
            parent = None
            for potential_parent in tree.get_nonterminals():
                if potential_parent != current:
                    for child in potential_parent.clades:
                        if child == current:
                            parent = potential_parent
                            break
                    if parent:
                        break
            
            if parent == node:
                is_descendant = True
                break
            elif parent is None:
                break
            current = parent
        
        if is_descendant:
            descendants.append(tip)
    
    # Clade similarity is defined only by descendant terminal sequences.
    for desc in descendants:
        if desc.name in tip_sequences:
            sequences.append(tip_sequences[desc.name])
    
    return sequences

def calculate_node_similarity(sequences):
    """
    计算一个节点包含的所有序列之间的最低相似性
    用于与阈值进行比较分析
    """
    if len(sequences) <= 1:
        return 100.0
    
    similarities = []
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            try:
                sim = calculate_sequence_similarity(sequences[i], sequences[j])
                similarities.append(sim)
            except ValueError:
                continue
    
    if not similarities:
        return 0.0
    
    return min(similarities)  # 使用最低相似性用于阈值比较

def find_optimal_nodes(tree, ancestral_sequences, tip_sequences, similarity_threshold=95.0):
    """
    寻找满足条件的节点：
    1. 节点内部序列的最低相似性 >= 阈值
    2. 该节点的parent节点内部序列相似性 < 阈值
    
    Args:
        tree: 系统发育树（已进行midpoint rooting和节点编号）
        ancestral_sequences (dict): 祖先序列字典
        tip_sequences (dict): 末端序列字典
        similarity_threshold (float): 相似性阈值
        
    Returns:
        list: 满足条件的节点列表，每个元素为(node, node_similarity, is_terminal_node)
    """
    root = tree.root
    optimal_nodes = []
    
    logging.info(f"开始分析树结构，根节点: {root.name if hasattr(root, 'name') else 'Root'}")
    logging.info(f"相似性阈值: {similarity_threshold}%")
    
    # 获取所有节点
    all_nodes = list(tree.get_nonterminals()) + list(tree.get_terminals())
    # 确保根节点在列表中
    if root not in all_nodes:
        all_nodes.insert(0, root)
    
    # 为每个节点计算相似性
    node_similarities = {}
    for node in all_nodes:
        node_seqs = get_node_sequences(tree, node, ancestral_sequences, tip_sequences)
        if node_seqs:
            similarity = calculate_node_similarity(node_seqs)
            node_similarities[node] = similarity
            logging.info(f"节点 {node.name}: {len(node_seqs)} 条序列, 相似性: {similarity:.2f}%")
        else:
            node_similarities[node] = 0.0
            logging.warning(f"节点 {node.name} 没有序列")
    
    # 为每个节点找到其parent节点
    node_parents = {}
    for node in all_nodes:
        if node == root:
            node_parents[node] = None
        else:
            # 找到包含该节点的最小内部节点作为parent
            parent = None
            for potential_parent in tree.get_nonterminals():
                if potential_parent != node:
                    for child in potential_parent.clades:
                        if child == node:
                            parent = potential_parent
                            break
                    if parent:
                        break
            node_parents[node] = parent
    
    # 检查每个节点是否满足条件
    for node in all_nodes:
        current_similarity = node_similarities[node]
        parent = node_parents[node]
        is_terminal = node.is_terminal()
        
        # 检查条件1: 当前节点相似性 >= 阈值
        condition1 = current_similarity >= similarity_threshold
        
        # 检查条件2: parent节点相似性 < 阈值（如果存在parent）
        condition2 = True
        if parent is not None:
            parent_similarity = node_similarities[parent]
            condition2 = parent_similarity < similarity_threshold
            logging.info(f"节点 {node.name}: 相似性 {current_similarity:.2f}% >= {similarity_threshold}%: {condition1}")
            logging.info(f"节点 {node.name}: parent {parent.name} 相似性 {parent_similarity:.2f}% < {similarity_threshold}%: {condition2}")
        else:
            logging.info(f"节点 {node.name}: 相似性 {current_similarity:.2f}% >= {similarity_threshold}%: {condition1}")
            logging.info(f"节点 {node.name}: 无parent节点")
        
        # 如果两个条件都满足，添加到结果中
        if condition1 and condition2:
            optimal_nodes.append((node, current_similarity, is_terminal))
            logging.info(f"✅ 找到满足条件的节点: {node.name}, 相似性: {current_similarity:.2f}%, 末端节点: {is_terminal}")
        else:
            if is_terminal and current_similarity < similarity_threshold:
                # 末端节点且不满足阈值，仍然保留
                optimal_nodes.append((node, current_similarity, is_terminal))
                logging.info(f"⚠️ 保留末端节点: {node.name}, 相似性: {current_similarity:.2f}% < {similarity_threshold}%, 但作为代表序列保留")
    
    # 如果没有找到任何节点，返回根节点
    if not optimal_nodes:
        root_seqs = get_node_sequences(tree, root, ancestral_sequences, tip_sequences)
        root_similarity = calculate_node_similarity(root_seqs)
        optimal_nodes.append((root, root_similarity, False))
        logging.info(f"⚠️ 未找到满足条件的节点，使用根节点: {root.name}, 相似性: {root_similarity:.2f}%")
    
    logging.info(f"总共找到 {len(optimal_nodes)} 个节点")
    return optimal_nodes

def main():
    parser = argparse.ArgumentParser(description="Enhanced Consensus Generation with Validation")
    parser.add_argument("--input_file", required=True, help="Input consensus FASTA file from Process 5")
    parser.add_argument("--output_file", required=True, help="Output enhanced consensus FASTA file")
    parser.add_argument("--similarity_threshold", type=float, default=95.0, help="Similarity threshold")
    parser.add_argument("--max_iterations", type=int, default=10, help="Maximum iterations")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--fast_mode", action='store_true', help="快速模式：减少bootstrap次数，使用更少的迭代")
    
    args = parser.parse_args()
    
    logging.info(f"开始enhanced consensus处理...")
    logging.info(f"输入文件: {args.input_file}")
    logging.info(f"输出文件: {args.output_file}")
    logging.info(f"相似性阈值: {args.similarity_threshold}%")
    logging.info(f"最大迭代次数: {args.max_iterations}")
    
    # 检查输入文件
    if not os.path.exists(args.input_file):
        logging.error(f"输入文件不存在: {args.input_file}")
        sys.exit(1)
    
    # 读取序列数量
    sequences = list(SeqIO.parse(args.input_file, "fasta"))
    if len(sequences) <= 1:
        logging.info("序列数量<=1，直接复制输入文件")
        with open(args.output_file, 'w') as f:
            SeqIO.write(sequences, f, "fasta")
        sys.exit(0)
    
    logging.info(f"输入文件包含 {len(sequences)} 条序列")
    
    # 创建临时目录进行迭代处理
    with tempfile.TemporaryDirectory() as temp_dir:
        input_stem = os.path.splitext(os.path.basename(args.input_file))[0]
        current_input = os.path.join(temp_dir, os.path.basename(args.input_file))
        SeqIO.write(
            make_record_ids_unique(list(SeqIO.parse(args.input_file, "fasta"))),
            current_input,
            "fasta",
        )
        current_iteration = 1
        
        while current_iteration <= args.max_iterations:
            logging.info(f"\n=== 第 {current_iteration} 轮处理 ===")
            
            # 步骤1: 计算当前文件的pairwise相似性
            max_similarity = calculate_maximum_pairwise_similarity(current_input)
            
            logging.info(f"当前最高相似性: {max_similarity:.4f}%")
            
            # 步骤2: 检查是否需要继续enhanced consensus
            if max_similarity < args.similarity_threshold:
                logging.info(f"相似性 {max_similarity:.4f}% < {args.similarity_threshold}%，满足条件，停止迭代")
                break
            
            logging.info(f"相似性 {max_similarity:.4f}% >= {args.similarity_threshold}%，需要继续enhanced consensus")
            
            # 步骤3: 运行generate_consensus.py的流程生成consensus序列
            current_output = os.path.join(
                temp_dir, f"{input_stem}.enhanced{current_iteration}.fasta"
            )
            
            try:
                # 调用generate_consensus.py的流程
                generate_consensus_sequences(
                    current_input, 
                    current_output, 
                    args.similarity_threshold,
                    args.threads,
                    current_iteration # 传递当前迭代次数
                )
                normalize_fasta_ids(current_output)
                logging.info(f"第{current_iteration}轮generate_consensus.py流程完成，输出文件: {current_output}")
                merged_pair_count = tree_guided_merge_similar_consensus(
                    current_output,
                    threshold=args.similarity_threshold,
                    threads=args.threads,
                    iteration=current_iteration,
                )
                if merged_pair_count:
                    normalize_fasta_ids(current_output)
                    logging.info(
                        "第%d轮树引导合并了 %d 对高相似consensus代表序列",
                        current_iteration,
                        merged_pair_count,
                    )
                
                # 检查输出文件是否存在
                if not os.path.exists(current_output):
                    logging.error(f"第{current_iteration}轮输出文件不存在: {current_output}")
                    break
                    
                # 检查输出文件的序列数量
                output_sequences = list(SeqIO.parse(current_output, "fasta"))
                logging.info(f"第{current_iteration}轮生成了 {len(output_sequences)} 条consensus序列")
                
                if len(output_sequences) <= 1:
                    current_input = current_output
                    current_iteration += 1
                    logging.info(f"第{current_iteration}轮只生成了 {len(output_sequences)} 条序列，停止迭代")
                    break
                    
            except Exception as e:
                logging.error(f"第{current_iteration}轮generate_consensus.py流程失败: {str(e)}")
                break
            
            # 更新输入文件为当前输出
            current_input = current_output
            current_iteration += 1
        
        # 复制最终结果到输出文件
        if current_iteration > args.max_iterations:
            logging.warning(f"达到最大迭代次数 {args.max_iterations}，使用最后结果")
        
        # current_input始终指向最后一次成功完成的迭代结果。
        final_output_file = current_input
        logging.info(f"使用最终迭代文件: {final_output_file}")
        
        with open(final_output_file, 'r') as src, open(args.output_file, 'w') as dst:
            dst.write(src.read())
        
        logging.info(f"Enhanced consensus处理完成，共进行 {current_iteration-1} 轮迭代")
        logging.info(f"最终结果保存到: {args.output_file}")

def generate_consensus_sequences(aligned_fasta, consensus_output, similarity_threshold=95.0, threads=8, iteration=1):
    """
    调用generate_consensus.py的流程生成consensus序列
    """
    logging.info(f"使用generate_consensus.py流程处理: {aligned_fasta}")
    
    # 读取比对文件
    try:
        alignment = list(SeqIO.parse(aligned_fasta, "fasta"))
        num_sequences = len(alignment)
    except Exception as e:
        logging.error(f"读取比对文件失败: {str(e)}")
        raise
    
    logging.info(f"检测到 {num_sequences} 条序列")
    
    # 从输入文件名提取基因型和聚类ID信息
    input_basename = os.path.basename(aligned_fasta)
    input_name = os.path.splitext(input_basename)[0]
    
    if '_consensus_' in input_name:
        parts = input_name.split('_consensus_')
        if len(parts) == 2:
            genotype = parts[0]
            cluster_id = parts[1]
        else:
            genotype = "Unknown"
            cluster_id = "Unknown"
    else:
        genotype = "Unknown"
        cluster_id = "Unknown"
    
    logging.info(f"从输入文件名解析: Genotype: {genotype}, Cluster ID: {cluster_id}")
    
    if num_sequences == 1:
        # 只有1条序列，直接保存为一致性序列
        record = alignment[0]
        
        # 从cluster序列中提取共有基因型
        common_genotype = extract_common_genotype_from_cluster(aligned_fasta)
        logging.info(f"从cluster序列中提取的共有基因型: {common_genotype}")
        
        # 从文件名中提取cluster序号
        cluster_number = extract_cluster_number_from_filename(aligned_fasta)
        logging.info(f"从文件名中提取的cluster序号: {cluster_number}")
        
        # 从文件名中提取enhanced迭代次数
        enhanced_iteration = extract_enhanced_iteration_from_filename(aligned_fasta)
        
        # 生成consensus序列ID：共有基因型_cluster序号_enhanced迭代次数_node_1（单条序列只有一个节点）
        if enhanced_iteration == 1:
            consensus_id = f"{common_genotype}_cluster_{cluster_number}_enhanced1_node_1"
        else:
            consensus_id = f"{common_genotype}_cluster_{cluster_number}_enhanced{iteration}_node_1"
        
        record.id = consensus_id
        record.description = f"Single sequence consensus for {common_genotype} cluster {cluster_number}"
        SeqIO.write(record, consensus_output, "fasta")
        logging.info(f"聚类 {common_genotype} cluster {cluster_number} 仅有1条序列，直接输出为一致性序列")
        
    elif num_sequences == 2:
        # 2条序列，使用简单方法生成一致性序列
        generate_simple_consensus(aligned_fasta, consensus_output, genotype, cluster_id, iteration)
        logging.info(f"聚类 {genotype} cluster {cluster_id} 包含2条序列，使用简单方法生成一致性序列")
        
    else:
        # 3条或更多序列，使用节点分析方法
        try:
            generate_consensus_with_node_analysis(
                aligned_fasta, consensus_output, genotype, cluster_id, 
                similarity_threshold, threads, iteration
            )
            logging.info(f"聚类 {genotype} cluster {cluster_id} 包含 {num_sequences} 条序列，使用节点分析方法生成一致性序列")
        except Exception as e:
            logging.error(f"节点分析方法失败，回退到简单方法: {str(e)}")
            generate_simple_consensus(aligned_fasta, consensus_output, genotype, cluster_id)
            logging.info(f"聚类 {genotype} cluster {cluster_id} 使用回退方法生成一致性序列")

def generate_simple_consensus(aligned_file, consensus_file, genotype, cluster_id, iteration=1):
    """生成简单的一致性序列（用于2条序列的情况）"""
    try:
        sequences = list(SeqIO.parse(aligned_file, "fasta"))
        if len(sequences) != 2:
            raise ValueError("简单一致性序列生成只适用于2条序列")
        
        # 从文件名中提取enhanced迭代次数
        enhanced_iteration = extract_enhanced_iteration_from_filename(aligned_file)
        
        # 计算两条序列的相似性
        seq1 = str(sequences[0].seq)
        seq2 = str(sequences[1].seq)
        similarity = calculate_sequence_similarity(seq1, seq2)
        logging.info(f"两条序列的相似性: {similarity:.2f}%")
        
        if similarity >= 95.0:
            # 只输出一条共识序列，ID包含enhanced迭代次数标记
            consensus_id = f"{genotype}_cluster_{cluster_id}_enhanced{iteration}_consensus"
            consensus_record = SeqRecord(
                Seq(seq1),
                id=consensus_id,
                description=""
            )
            SeqIO.write(consensus_record, consensus_file, "fasta")
            logging.info(f"两条序列相似性 >= 95%，只输出一条共识序列: {consensus_id}")
        else:
            # 如果相似性<95%，两条序列都保留
            logging.info(f"两条序列相似性 < 95%，两条序列都保留")
            with open(consensus_file, "w") as f:
                consensus_id_1 = f"{genotype}_cluster_{cluster_id}_enhanced{iteration}_node_1"
                consensus_id_2 = f"{genotype}_cluster_{cluster_id}_enhanced{iteration}_node_2"
                f.write(f">{consensus_id_1}\n{seq1}\n")
                f.write(f">{consensus_id_2}\n{seq2}\n")
            logging.info(f"两条序列已保存至: {consensus_file}")
    except Exception as e:
        logging.error(f"生成简单一致性序列失败: {str(e)}")
        raise

def generate_consensus_with_node_analysis(aligned_file, consensus_file, genotype, cluster_id, similarity_threshold=95.0, threads=8, iteration=1):
    """使用节点分析方法生成一致性序列"""
    logging.info(f"使用节点分析生成一致性序列: {aligned_file}")
    
    # 从cluster序列中提取共有基因型
    common_genotype = extract_common_genotype_from_cluster(aligned_file)
    logging.info(f"从cluster序列中提取的共有基因型: {common_genotype}")
    
    # 从文件名中提取cluster序号
    cluster_number = extract_cluster_number_from_filename(aligned_file)
    logging.info(f"从文件名中提取的cluster序号: {cluster_number}")
    
    # 从文件名中提取enhanced迭代次数
    enhanced_iteration = extract_enhanced_iteration_from_filename(aligned_file)
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as temp_dir:
        # 步骤1: 使用IQ-TREE构建系统发育树
        tree_prefix = os.path.join(temp_dir, "tree")
        run_iqtree_enhanced(aligned_file, tree_prefix, threads)
        
        tree_file = f"{tree_prefix}.treefile"
        ancestral_file = f"{tree_prefix}.state"
        
        # 步骤2: 对树进行midpoint rooting
        rooted_tree_file = os.path.join(temp_dir, "tree_rooted.treefile")
        midpoint_root_tree_manual(tree_file, rooted_tree_file)
        
        # 步骤3: 为树的节点分配编号
        tree = Phylo.read(rooted_tree_file, "newick")
        assign_node_numbers(tree)
        
        # 步骤4: 解析祖先序列和末端序列
        ancestral_sequences = parse_ancestral_sequences(ancestral_file)
        tip_sequences = parse_tip_sequences(aligned_file)
        
        # 步骤5: 寻找最优节点
        optimal_nodes = find_optimal_nodes(tree, ancestral_sequences, tip_sequences, similarity_threshold)
        
        logging.info(f"找到 {len(optimal_nodes)} 个满足条件的节点")
        
        # 步骤6: 为每个节点生成consensus序列
        consensus_records = []
        for i, (node, node_similarity, is_terminal) in enumerate(optimal_nodes):
            logging.info(f"处理节点 {i+1}: {node.name}, 节点最小相似性: {node_similarity:.2f}%, 是否为末端节点: {is_terminal}")
            
            # 获取节点的基因型信息
            node_genotypes = get_node_genotypes(tree, node, ancestral_sequences, tip_sequences)
            node_genotype_name = get_node_genotype_name(node_genotypes)
            logging.info(f"节点 {node.name} 包含的基因型: {node_genotypes}, 节点基因型名称: {node_genotype_name}")
            
            # 获取节点的序列
            node_sequences = get_node_sequences(tree, node, ancestral_sequences, tip_sequences)
            
            if not node_sequences:
                logging.warning(f"节点 {node.name} 没有序列，跳过")
                continue
            
            # 生成consensus序列
            consensus_seq = generate_simple_consensus_from_sequences(node_sequences)
            
            # 根据节点类型和迭代次数生成序列ID
            if is_terminal and len(node_sequences) == 1:
                # 末端节点且只有1条序列，保持原始名称
                consensus_id = node.name
            else:
                # 内部节点或多条序列，添加enhanced标识
                # 使用从文件名中提取的enhanced迭代次数，而不是传入的iteration参数
                consensus_id = f"{node_genotype_name}_cluster_{cluster_number}_enhanced{iteration}_node_{i+1}"
            
            # 创建consensus记录
            consensus_record = SeqRecord(
                Seq(consensus_seq),
                id=consensus_id,
                description=""
            )
            
            consensus_records.append(consensus_record)
            logging.info(f"生成consensus序列: {consensus_id}")
        
        # 保存所有consensus序列
        SeqIO.write(make_record_ids_unique(consensus_records), consensus_file, "fasta")
        logging.info(f"共生成 {len(consensus_records)} 个consensus序列，已保存至: {consensus_file}")

if __name__ == "__main__":
    main() 
