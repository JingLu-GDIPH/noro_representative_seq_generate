#!/usr/bin/env python3
"""
Generate consensus sequences using phylogenetic tree reconstruction with node analysis.

This script processes aligned FASTA files and generates consensus sequences
by analyzing all nodes in the phylogenetic tree and selecting the optimal node
that meets similarity threshold criteria.
"""

import os
import sys
import argparse
import csv
import subprocess
import tempfile
import logging
import re
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
import numpy as np

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

VALID_BASES = set("ACGT")
VALID_BASE_BYTES = np.frombuffer(b"ACGT", dtype=np.uint8)
RDRP_GENOTYPE_PATTERN = re.compile(r"^G(?:I|II|IX)\.P[A-Za-z0-9]+$")
VP1_GENOTYPE_PATTERN = re.compile(r"^G(?:I|II|IX)\.[A-Za-z0-9]+$")
GENOTYPE_PAIR_PATTERN = re.compile(
    r"(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)"
)


def extract_genotype_pair_label(text):
    """Extract a combined RdRp_VP1 genotype label when present."""
    match = GENOTYPE_PAIR_PATTERN.search(str(text))
    if match:
        return f"{match.group(1)}_{match.group(2)}"
    return None


def extract_genotype_pair_from_id(sequence_id):
    """Extract RdRp_VP1 from the first two FASTA ID fields."""
    parts = str(sequence_id).split("_")
    if len(parts) >= 2 and RDRP_GENOTYPE_PATTERN.match(parts[0]) and VP1_GENOTYPE_PATTERN.match(parts[1]):
        return f"{parts[0]}_{parts[1]}"
    return extract_genotype_pair_label(sequence_id)

def calculate_sequence_similarity(seq1, seq2):
    """
    计算两条序列的相似性。

    只有当两个位置都为 A/C/G/T 时才计入分母；N 和 gap 均不作为有效覆盖。
    
    Args:
        seq1, seq2 (str): 两条序列
        
    Returns:
        float: 相似性百分比 (0-100)
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length")
    
    if len(seq1) == 0:
        return 0.0
    
    matches = 0
    total_positions = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a not in VALID_BASES or b not in VALID_BASES:
            continue
        total_positions += 1
        if a == b:
            matches += 1
    
    if total_positions == 0:
        return 0.0
    
    return (matches / total_positions) * 100

def calculate_node_similarity(sequences):
    """
    计算一个节点包含的所有序列之间的最小相似性
    
    Args:
        sequences (list): 序列列表
        
    Returns:
        float: 最小相似性百分比
    """
    if len(sequences) <= 1:
        return 100.0

    sequence_lengths = {len(sequence) for sequence in sequences}
    if len(sequence_lengths) != 1:
        return 0.0

    seq_length = sequence_lengths.pop()
    if seq_length == 0:
        return 0.0

    encoded = np.frombuffer(
        "".join(sequence.upper() for sequence in sequences).encode("ascii"),
        dtype=np.uint8
    ).reshape(len(sequences), seq_length)
    min_similarity = 100.0
    valid_pairs = 0

    # 每次将一条序列与其后的所有序列同时比较，避免Python逐碱基循环。
    for i in range(len(sequences) - 1):
        reference = encoded[i]
        comparisons = encoded[i + 1:]
        valid = np.isin(comparisons, VALID_BASE_BYTES) & np.isin(reference, VALID_BASE_BYTES)
        totals = valid.sum(axis=1)
        matches = ((comparisons == reference) & valid).sum(axis=1)
        nonzero = totals > 0

        if np.any(nonzero):
            similarities = matches[nonzero] / totals[nonzero] * 100.0
            min_similarity = min(min_similarity, float(similarities.min()))
            valid_pairs += int(nonzero.sum())

    return min_similarity if valid_pairs else 0.0

def get_node_sequences(tree, node, ancestral_sequences, tip_sequences):
    """
    获取指定节点包含的所有序列（包括祖先序列和后代序列）
    
    Args:
        tree: 系统发育树
        node: 目标节点
        ancestral_sequences (dict): 祖先序列字典
        tip_sequences (dict): 末端序列字典
        
    Returns:
        list: 该节点包含的所有序列
    """
    sequences = []
    
    # 如果是末端节点，直接返回该节点的序列
    if node.is_terminal():
        if node.name in tip_sequences:
            sequences.append(tip_sequences[node.name])
        return sequences
    
    # Clade similarity is defined only by descendant terminal sequences.
    for desc in node.get_terminals():
        if desc.name in tip_sequences:
            sequences.append(tip_sequences[desc.name])
    
    return sequences

def has_pairwise_below_threshold(sequences, threshold):
    """
    判断序列集合中是否存在任意一对pairwise相似性低于阈值
    Args:
        sequences (list): 序列列表
        threshold (float): 阈值
    Returns:
        bool: 是否存在一对pairwise相似性低于阈值
    """
    if len(sequences) <= 1:
        return False
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            try:
                sim = calculate_sequence_similarity(sequences[i], sequences[j])
                if sim < threshold:
                    return True
            except ValueError:
                continue
    return False

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

def run_iqtree(aligned_file, output_prefix, threads=8):
    """使用IQ-TREE构建系统发育树并生成祖先序列"""
    logging.info(f"使用IQ-TREE构建系统发育树和祖先序列: {aligned_file}")
    
    cmd = [
        'iqtree',
        '-s', aligned_file,
        '-m', 'JC',
        '-T', str(threads),
        '-redo',
        '-asr',  # 启用祖先序列重建
        '-pre', output_prefix
    ]
    
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"IQ-TREE建树和祖先序列重建完成，结果保存在 {output_prefix}.treefile")
    except subprocess.CalledProcessError as e:
        logging.error(f"IQ-TREE运行失败: {str(e)}")
        raise

def midpoint_root_tree(tree_file, output_file):
    """对树进行midpoint rooting"""
    logging.info(f"对树进行midpoint rooting: {tree_file}")
    
    # 使用IQ-TREE的midpoint rooting功能
    cmd = [
        'iqtree',
        '-t', tree_file,
        '-midpoint',
        '-pre', output_file.replace('.treefile', '')
    ]
    
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"Midpoint rooting完成，结果保存在: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Midpoint rooting失败: {str(e)}")
        # 如果失败，直接复制原文件
        import shutil
        shutil.copy2(tree_file, output_file)
        logging.warning(f"使用原始树文件: {output_file}")

def midpoint_root_tree_manual(tree_file, output_file):
    """手动实现midpoint rooting"""
    logging.info(f"手动进行midpoint rooting: {tree_file}")
    
    try:
        from Bio import Phylo
        
        # 读取树
        tree = Phylo.read(tree_file, "newick")
        
        # 计算树的总长度
        total_length = tree.total_branch_length()
        
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
            Phylo.write(tree, output_file, "newick")
            logging.info(f"Midpoint rooting完成，结果保存在: {output_file}")
        else:
            # 如果无法找到最长路径，直接复制原文件
            import shutil
            shutil.copy2(tree_file, output_file)
            logging.warning(f"无法找到最长路径，使用原始树文件: {output_file}")
            
    except Exception as e:
        logging.error(f"手动midpoint rooting失败: {str(e)}")
        # 如果失败，直接复制原文件
        import shutil
        shutil.copy2(tree_file, output_file)
        logging.warning(f"使用原始树文件: {output_file}")

def assign_node_numbers(tree):
    """为树的节点分配编号"""
    logging.info("为树的节点分配编号")
    
    node_counter = 1
    
    # 为所有内部节点分配编号
    for node in tree.get_nonterminals():
        if not hasattr(node, 'name') or not node.name:
            node.name = f"Node{node_counter}"
            node_counter += 1
        elif not node.name.startswith('Node'):
            # 如果节点有名称但不是Node格式，保留原名称
            pass
    
    logging.info(f"为 {node_counter-1} 个内部节点分配了编号")
    return tree

def parse_ancestral_sequences(ancestral_file):
    """
    解析祖先序列文件
    
    Args:
        ancestral_file (str): 祖先序列文件路径
        
    Returns:
        dict: 节点名称到序列的映射
    """
    if not os.path.exists(ancestral_file):
        logging.warning(f"祖先序列文件不存在: {ancestral_file}")
        return {}

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
        logging.info(f"成功解析 {len(ancestral_sequences)} 个祖先序列")
        return ancestral_sequences
    except Exception as e:
        logging.error(f"解析祖先序列文件失败: {str(e)}")
        return {}

def parse_tip_sequences(aligned_file):
    """
    解析末端序列文件
    
    Args:
        aligned_file (str): 比对后的FASTA文件路径
        
    Returns:
        dict: 序列名称到序列的映射
    """
    tip_sequences = {}
    
    try:
        for record in SeqIO.parse(aligned_file, "fasta"):
            tip_sequences[record.id] = str(record.seq)
        
        logging.info(f"成功解析 {len(tip_sequences)} 个末端序列")
    except Exception as e:
        logging.error(f"解析末端序列文件失败: {str(e)}")
    
    return tip_sequences

def extract_common_genotype_from_cluster(aligned_file):
    """
    从cluster的序列中提取共有的基因型
    如果所有序列都是同一基因型，返回该基因型；否则返回"mixed"
    
    Args:
        aligned_file (str): 比对后的FASTA文件路径
        
    Returns:
        str: 共有基因型或"mixed"
    """
    genotypes = set()
    
    try:
        filename_genotype = extract_genotype_pair_label(os.path.basename(aligned_file))
        if filename_genotype:
            return filename_genotype

        for record in SeqIO.parse(aligned_file, "fasta"):
            genotype = extract_genotype_pair_from_id(record.id)
            if genotype:
                genotypes.add(genotype)
    except Exception as e:
        logging.error(f"解析序列基因型失败: {str(e)}")
        return "unknown"
    
    if len(genotypes) == 1:
        return list(genotypes)[0]
    else:
        return "mixed"

def extract_cluster_number_from_filename(filename):
    """
    从文件名中提取cluster序号
    
    Args:
        filename (str): 文件名
        
    Returns:
        str: cluster序号
    """
    # 尝试从文件名中提取cluster信息
    basename = os.path.basename(filename)
    name_without_ext = os.path.splitext(basename)[0]
    
    # 查找cluster_数字的模式
    import re
    cluster_match = re.search(r'cluster_(\d+)', name_without_ext)
    if cluster_match:
        return cluster_match.group(1)
    
    # 如果没有找到，尝试其他模式
    parts = name_without_ext.split('_')
    for i, part in enumerate(parts):
        if part == 'cluster' and i + 1 < len(parts):
            return parts[i + 1]
    
    return "unknown"

def extract_enhanced_iteration_from_filename(filename):
    """
    从文件名中提取enhanced consensus的迭代次数
    
    Args:
        filename (str): 文件名
        
    Returns:
        int: 迭代次数，如果是第一次则为1
    """
    basename = os.path.basename(filename)
    name_without_ext = os.path.splitext(basename)[0]
    
    # 检查是否包含enhanced标识
    if 'enhanced' in name_without_ext.lower():
        # 尝试提取enhanced后的数字
        import re
        enhanced_match = re.search(r'enhanced(\d+)', name_without_ext.lower())
        if enhanced_match:
            return int(enhanced_match.group(1))
        else:
            # 如果只有enhanced但没有数字，认为是第1次
            return 1
    
    # 如果没有enhanced标识，认为是第1次
    return 1

def generate_consensus_with_node_analysis(aligned_file, consensus_file, genotype, cluster_id, 
                                        similarity_threshold=95.0, threads=8):
    """
    使用节点分析生成一致性序列
    
    Args:
        aligned_file (str): 比对后的FASTA文件路径
        consensus_file (str): 输出一致性序列文件路径
        genotype (str): 基因型（可能被覆盖）
        cluster_id (str): 聚类ID
        similarity_threshold (float): 相似性阈值
        threads (int): 线程数
    """
    logging.info(f"使用节点分析生成一致性序列: {aligned_file}")
    
    # 从cluster序列中提取共有基因型
    common_genotype = extract_common_genotype_from_cluster(aligned_file)
    logging.info(f"从cluster序列中提取的共有基因型: {common_genotype}")
    
    # 从文件名中提取cluster序号
    cluster_number = extract_cluster_number_from_filename(aligned_file)
    logging.info(f"从文件名中提取的cluster序号: {cluster_number}")
    
    # 从文件名中提取enhanced迭代次数
    enhanced_iteration = extract_enhanced_iteration_from_filename(aligned_file)
    logging.info(f"从文件名中提取的enhanced迭代次数: {enhanced_iteration}")
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as tmpdir:
        # 运行IQ-TREE构建系统发育树
        tree_prefix = os.path.join(tmpdir, "tree")
        run_iqtree(aligned_file, tree_prefix, threads)
        
        # 读取原始系统发育树
        original_tree_file = f"{tree_prefix}.treefile"
        if not os.path.exists(original_tree_file):
            raise FileNotFoundError(f"树文件不存在: {original_tree_file}")
        
        # 进行midpoint rooting
        rooted_tree_file = f"{tree_prefix}_rooted.treefile"
        midpoint_root_tree_manual(original_tree_file, rooted_tree_file)
        
        # 读取rooted树
        tree = Phylo.read(rooted_tree_file, "newick")
        
        # 为树分配节点编号
        tree = assign_node_numbers(tree)
        
        # 解析祖先序列
        ancestral_file = f"{tree_prefix}.state"
        ancestral_sequences = parse_ancestral_sequences(ancestral_file)
        
        # 解析末端序列
        tip_sequences = parse_tip_sequences(aligned_file)
        
        # 找到所有满足条件的节点
        optimal_nodes = find_optimal_nodes(tree, ancestral_sequences, tip_sequences, similarity_threshold)
        
        logging.info(f"找到 {len(optimal_nodes)} 个满足条件的节点")
        
        # 为每个满足条件的节点生成consensus序列
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
                    # 使用简单方法生成一致性序列
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
            # 对于末端节点，我们需要检查其父节点是否满足阈值
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
                    # 如果无法获取原始名称，使用默认格式（Process 5不包含enhanced1）
                    consensus_id = f"{node_genotype_name}_cluster_{cluster_number}_node_{node_num}"
            else:
                # 对于满足条件的节点或内部节点，使用新的命名格式（Process 5不包含enhanced1）
                consensus_id = f"{node_genotype_name}_cluster_{cluster_number}_node_{node_num}"
            
            consensus_sequences.append((consensus_id, optimal_sequence))
            logging.info(f"生成consensus序列: {consensus_id}")
        
        # 保存所有consensus序列
        with open(consensus_file, "w") as f:
            for consensus_id, sequence in consensus_sequences:
                f.write(f">{consensus_id}\n{sequence}\n")
        
        logging.info(f"共生成 {len(consensus_sequences)} 个consensus序列，已保存至: {consensus_file}")

def generate_simple_consensus_from_sequences(sequences):
    """
    从序列列表生成简单的一致性序列
    
    Args:
        sequences (list): 序列列表
        
    Returns:
        str: 一致性序列
    """
    if not sequences:
        raise ValueError("序列列表为空")
    
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
    
    # 计算每个位置的最常见碱基
    consensus = []
    seq_length = len(sequences[0])
    
    for pos in range(seq_length):
        bases = [
            seq[pos].upper()
            for seq in sequences
            if pos < len(seq) and seq[pos].upper() in VALID_BASES
        ]
        
        if not bases:
            consensus.append('-')
        else:
            # 计算最常见的碱基
            from collections import Counter
            base_counts = Counter(bases)
            most_common = base_counts.most_common(1)[0][0]
            consensus.append(most_common)
    
    return ''.join(consensus)

def generate_simple_consensus(aligned_file, consensus_file, genotype, cluster_id):
    """使用简单方法生成一致性序列（适用于2条序列）"""
    logging.info(f"使用简单方法生成一致性序列: {aligned_file}")
    
    # 从cluster序列中提取共有基因型
    common_genotype = extract_common_genotype_from_cluster(aligned_file)
    logging.info(f"从cluster序列中提取的共有基因型: {common_genotype}")
    
    # 从文件名中提取cluster序号
    cluster_number = extract_cluster_number_from_filename(aligned_file)
    logging.info(f"从文件名中提取的cluster序号: {cluster_number}")
    
    # 从文件名中提取enhanced迭代次数
    enhanced_iteration = extract_enhanced_iteration_from_filename(aligned_file)
    logging.info(f"从文件名中提取的enhanced迭代次数: {enhanced_iteration}")
    
    alignment = AlignIO.read(aligned_file, "fasta")
    
    # 检查是否有两条序列
    if len(alignment) == 2:
        seq1 = str(alignment[0].seq)
        seq2 = str(alignment[1].seq)
        
        # 计算两条序列的相似性
        similarity = calculate_sequence_similarity(seq1, seq2)
        logging.info(f"两条序列的相似性: {similarity:.2f}%")
        
        # 如果相似性>=95%，选择第一条序列作为代表
        if similarity >= 95.0:
            logging.info(f"两条序列相似性 >= 95%，选择第一条序列作为代表")
            selected_sequence = seq1
            consensus_id = f"{common_genotype}_cluster_{cluster_number}_node_1"
            
            with open(consensus_file, "w") as f:
                f.write(f">{consensus_id}\n{selected_sequence}\n")
            logging.info(f"单条序列已保存至: {consensus_file}")
            return
        else:
            # 如果相似性<95%，两条序列都保留
            logging.info(f"两条序列相似性 < 95%，两条序列都保留")
            with open(consensus_file, "w") as f:
                # 生成consensus序列ID：共有基因型_cluster序号_node_1
                consensus_id_1 = f"{common_genotype}_cluster_{cluster_number}_node_1"
                consensus_id_2 = f"{common_genotype}_cluster_{cluster_number}_node_2"
                f.write(f">{consensus_id_1}\n{seq1}\n")
                f.write(f">{consensus_id_2}\n{seq2}\n")
            logging.info(f"两条序列已保存至: {consensus_file}")
            return
    else:
        # 其他情况使用原来的方法
        selected_sequence = generate_simple_consensus_from_sequences(
            [str(record.seq) for record in alignment]
        )
        consensus_id = f"{common_genotype}_cluster_{cluster_number}_node_1"
        
        with open(consensus_file, "w") as f:
            f.write(f">{consensus_id}\n{selected_sequence}\n")
        logging.info(f"简单一致性序列已保存至: {consensus_file}")

def generate_consensus_sequences(aligned_fasta, consensus_output, similarity_threshold=95.0, threads=8):
    """
    根据序列数量选择不同的方法生成一致性序列
    
    Args:
        aligned_fasta (str): 比对后的FASTA文件路径
        consensus_output (str): 输出一致性序列文件路径
        similarity_threshold (float): 相似性阈值
        threads (int): 使用的线程数
    """
    logging.info(f"开始处理比对文件: {aligned_fasta}")
    
    # 读取比对后的序列
    alignment = AlignIO.read(aligned_fasta, "fasta")
    num_sequences = len(alignment)
    
    logging.info(f"检测到 {num_sequences} 条序列")
    
    # 从输出文件名提取基因型和聚类ID信息
    output_basename = os.path.basename(consensus_output)
    output_name = os.path.splitext(output_basename)[0]
    
    if '_consensus_' in output_name:
        parts = output_name.split('_consensus_')
        if len(parts) == 2:
            genotype = parts[0]
            cluster_id = parts[1]
        else:
            genotype = "Unknown"
            cluster_id = "Unknown"
    else:
        genotype = "Unknown"
        cluster_id = "Unknown"
    
    logging.info(f"从输出文件名解析: Genotype: {genotype}, Cluster ID: {cluster_id}")
    
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
        
        # 生成consensus序列ID：共有基因型_cluster序号_node_1（Process 5不包含enhanced1）
        consensus_id = f"{common_genotype}_cluster_{cluster_number}_node_1"
        
        record.id = consensus_id
        record.description = f"Single sequence consensus for {common_genotype} cluster {cluster_number}"
        SeqIO.write(record, consensus_output, "fasta")
        logging.info(f"聚类 {common_genotype} cluster {cluster_number} 仅有1条序列，直接输出为一致性序列")
        
    elif num_sequences == 2:
        # 2条序列，使用简单方法生成一致性序列
        generate_simple_consensus(aligned_fasta, consensus_output, genotype, cluster_id)
        logging.info(f"聚类 {genotype} cluster {cluster_id} 包含2条序列，使用简单方法生成一致性序列")
        
    else:
        # 3条或更多序列，使用节点分析方法
        try:
            generate_consensus_with_node_analysis(
                aligned_fasta, consensus_output, genotype, cluster_id, 
                similarity_threshold, threads
            )
            logging.info(f"聚类 {genotype} cluster {cluster_id} 包含 {num_sequences} 条序列，使用节点分析方法生成一致性序列")
        except Exception as e:
            logging.error(f"节点分析方法失败，回退到简单方法: {str(e)}")
            generate_simple_consensus(aligned_fasta, consensus_output, genotype, cluster_id)
            logging.info(f"聚类 {genotype} cluster {cluster_id} 使用回退方法生成一致性序列")

def get_node_genotypes(tree, node, ancestral_sequences, tip_sequences):
    """
    获取指定节点包含的所有序列对应的基因型信息
    
    Args:
        tree: 系统发育树
        node: 目标节点
        ancestral_sequences (dict): 祖先序列字典
        tip_sequences (dict): 末端序列字典
        
    Returns:
        dict: 基因型到序列数量的映射
    """
    genotypes = {}
    
    # 如果是末端节点，直接从节点名称获取基因型
    if node.is_terminal():
        if node.name in tip_sequences:
            genotype = extract_genotype_pair_from_id(node.name) or node.name.split('_')[0]
            genotypes[genotype] = 1
        return genotypes
    
    # 统计所有后代序列的基因型
    for desc in node.get_terminals():
        if desc.name in tip_sequences:
            genotype = extract_genotype_pair_from_id(desc.name) or desc.name.split('_')[0]
            genotypes[genotype] = genotypes.get(genotype, 0) + 1
    
    return genotypes

def get_node_genotype_name(genotypes):
    """
    根据基因型统计生成节点基因型名称
    
    Args:
        genotypes (dict): 基因型到序列数量的映射
        
    Returns:
        str: 节点基因型名称
    """
    if not genotypes:
        return "unknown"
    
    # 如果只有一个基因型，直接返回
    if len(genotypes) == 1:
        return list(genotypes.keys())[0]
    
    # 如果有多个基因型，按数量排序并组合
    sorted_genotypes = sorted(genotypes.items(), key=lambda x: x[1], reverse=True)
    
    # 如果前两个基因型的数量相同且远大于其他基因型，使用前两个
    if len(sorted_genotypes) >= 2:
        top1_count = sorted_genotypes[0][1]
        top2_count = sorted_genotypes[1][1]
        total_count = sum(genotypes.values())
        
        # 如果前两个基因型占总数的80%以上，使用这两个基因型
        if (top1_count + top2_count) / total_count >= 0.8:
            return f"{sorted_genotypes[0][0]}_{sorted_genotypes[1][0]}"
    
    # 否则使用数量最多的基因型
    return sorted_genotypes[0][0]

def main():
    parser = argparse.ArgumentParser(description="Generate consensus sequences using phylogenetic node analysis.")
    parser.add_argument("--aligned_fasta", type=str, required=True, 
                       help="Path to the aligned FASTA file")
    parser.add_argument("--consensus_output", type=str, required=True, 
                       help="Path for the output consensus sequence file")
    parser.add_argument("--similarity_threshold", type=float, default=95.0,
                       help="Similarity threshold for node selection (default: 95.0)")
    parser.add_argument("--threads", type=int, default=8, 
                       help="Number of threads to use")
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.aligned_fasta):
        logging.error(f"输入文件不存在: {args.aligned_fasta}")
        sys.exit(1)
    
    # 检查必要的工具是否可用
    try:
        subprocess.run(['iqtree', '--version'], capture_output=True, check=True)
        logging.info("IQ-TREE 可用")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.warning("IQ-TREE 不可用，将使用简单方法")
    
    # 生成一致性序列
    try:
        generate_consensus_sequences(
            args.aligned_fasta, args.consensus_output, 
            args.similarity_threshold, args.threads
        )
        logging.info(f"一致性序列生成完成: {args.consensus_output}")
    except Exception as e:
        logging.error(f"生成一致性序列时发生错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 
