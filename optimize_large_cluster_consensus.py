#!/usr/bin/env python3
"""
优化大cluster的consensus生成
对于包含大量序列(>1000条)的cluster,采用分层采样策略:
1. 如果序列数<=500,使用原始方法
2. 如果序列数>500且<=1000,使用中等采样
3. 如果序列数>1000,使用分层聚类+采样策略
"""

import os
import sys
import argparse
import logging
import tempfile
import subprocess
import random
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_sequence_similarity(seq1, seq2):
    """计算两条序列的相似性"""
    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
    
    if len(seq1) == 0:
        return 0.0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-' and b != '-')
    total = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    
    if total == 0:
        return 0.0
    
    return (matches / total) * 100

def representative_sampling(sequences, sample_size, method='diverse'):
    """
    代表性采样
    
    Args:
        sequences: 序列列表
        sample_size: 采样大小
        method: 采样方法 ('random', 'diverse', 'stratified')
    
    Returns:
        采样后的序列列表
    """
    if len(sequences) <= sample_size:
        return sequences
    
    if method == 'random':
        # 随机采样
        return random.sample(sequences, sample_size)
    
    elif method == 'diverse':
        # 多样性采样:选择相互之间差异较大的序列
        selected = []
        remaining = sequences.copy()
        
        # 先随机选择第一条
        first = random.choice(remaining)
        selected.append(first)
        remaining.remove(first)
        
        # 迭代选择与已选序列差异最大的序列
        while len(selected) < sample_size and remaining:
            max_min_dist = -1
            best_seq = None
            
            # 从剩余序列中随机采样一部分进行评估(加速)
            sample_candidates = random.sample(remaining, min(100, len(remaining)))
            
            for candidate in sample_candidates:
                # 计算与所有已选序列的最小距离
                min_dist = min(
                    100 - calculate_sequence_similarity(str(candidate.seq), str(sel.seq))
                    for sel in selected
                )
                
                if min_dist > max_min_dist:
                    max_min_dist = min_dist
                    best_seq = candidate
            
            if best_seq:
                selected.append(best_seq)
                remaining.remove(best_seq)
            else:
                break
        
        # 如果还不够,随机补充
        if len(selected) < sample_size:
            needed = sample_size - len(selected)
            selected.extend(random.sample(remaining, min(needed, len(remaining))))
        
        return selected
    
    elif method == 'stratified':
        # 分层采样:按基因型分层
        genotype_groups = {}
        for seq in sequences:
            genotype = seq.id.split('_')[0]
            if genotype not in genotype_groups:
                genotype_groups[genotype] = []
            genotype_groups[genotype].append(seq)
        
        # 按比例从各基因型采样
        selected = []
        total_seqs = len(sequences)
        
        for genotype, seqs in genotype_groups.items():
            group_size = len(seqs)
            group_sample_size = max(1, int(sample_size * group_size / total_seqs))
            
            if group_size <= group_sample_size:
                selected.extend(seqs)
            else:
                selected.extend(random.sample(seqs, group_sample_size))
        
        # 如果采样不够,随机补充
        if len(selected) < sample_size:
            remaining = [s for s in sequences if s not in selected]
            needed = min(sample_size - len(selected), len(remaining))
            selected.extend(random.sample(remaining, needed))
        
        # 如果采样过多,随机删减
        if len(selected) > sample_size:
            selected = random.sample(selected, sample_size)
        
        return selected

def hierarchical_clustering_consensus(aligned_file, output_file, max_cluster_size=500, threads=4):
    """
    使用分层聚类生成consensus
    
    1. 先用vsearch对序列进行快速聚类
    2. 对每个子cluster生成consensus
    3. 合并所有sub-consensus
    """
    logging.info("使用分层聚类策略生成consensus")
    
    try:
        sequences = list(SeqIO.parse(aligned_file, "fasta"))
        num_sequences = len(sequences)
        
        logging.info(f"序列数: {num_sequences}")
        
        if num_sequences <= max_cluster_size:
            logging.info(f"序列数 <= {max_cluster_size}, 使用标准方法")
            return False
        
        # 创建临时目录
        with tempfile.TemporaryDirectory() as tmpdir:
            # 使用vsearch进行快速聚类
            temp_input = os.path.join(tmpdir, "input.fasta")
            temp_clusters = os.path.join(tmpdir, "clusters.uc")
            
            # 移除gap后保存(vsearch不接受gap)
            ungapped_sequences = []
            for seq in sequences:
                ungapped_seq = str(seq.seq).replace('-', '')
                if ungapped_seq:  # 确保不是空序列
                    ungapped_sequences.append(
                        SeqRecord(Seq(ungapped_seq), id=seq.id, description='')
                    )
            
            SeqIO.write(ungapped_sequences, temp_input, "fasta")
            
            # 计算合适的聚类阈值(较高的相似度以保持cluster内的一致性)
            identity_threshold = 0.97
            
            cmd = [
                'vsearch',
                '--cluster_fast', temp_input,
                '--id', str(identity_threshold),
                '--uc', temp_clusters,
                '--threads', str(threads),
                '--quiet'
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True, timeout=600)
            except subprocess.TimeoutExpired:
                logging.warning("vsearch聚类超时,回退到标准方法")
                return False
            except subprocess.CalledProcessError:
                logging.warning("vsearch聚类失败,回退到标准方法")
                return False
            
            # 解析聚类结果
            sub_clusters = {}
            with open(temp_clusters) as f:
                for line in f:
                    if line.startswith('H') or line.startswith('S'):
                        parts = line.strip().split('\t')
                        cluster_id = parts[1]
                        seq_id = parts[8]
                        
                        if cluster_id not in sub_clusters:
                            sub_clusters[cluster_id] = []
                        sub_clusters[cluster_id].append(seq_id)
            
            logging.info(f"分层聚类产生 {len(sub_clusters)} 个子clusters")
            
            # 为每个子cluster生成consensus(使用原始aligned序列)
            seq_dict = {seq.id: seq for seq in sequences}
            consensus_sequences = []
            
            for i, (cluster_id, seq_ids) in enumerate(sub_clusters.items(), 1):
                cluster_sequences = [seq_dict[sid] for sid in seq_ids if sid in seq_dict]
                
                if not cluster_sequences:
                    continue
                
                logging.info(f"处理子cluster {i}/{len(sub_clusters)}: {len(cluster_sequences)} 条序列")
                
                # 对子cluster生成consensus
                consensus_seq = generate_simple_consensus(cluster_sequences)
                
                consensus_record = SeqRecord(
                    Seq(consensus_seq),
                    id=f"subcluster_{i}_consensus",
                    description=f"Consensus of {len(cluster_sequences)} sequences"
                )
                consensus_sequences.append(consensus_record)
            
            # 保存所有sub-consensus
            if consensus_sequences:
                SeqIO.write(consensus_sequences, output_file, "fasta")
                logging.info(f"生成 {len(consensus_sequences)} 个sub-consensus序列")
                return True
            else:
                logging.warning("未生成任何sub-consensus")
                return False
    
    except Exception as e:
        logging.error(f"分层聚类失败: {str(e)}")
        return False

def generate_simple_consensus(sequences):
    """
    从序列列表生成简单的consensus序列
    使用多数原则
    """
    if not sequences:
        return ""
    
    if len(sequences) == 1:
        return str(sequences[0].seq)
    
    # 获取序列长度
    seq_length = len(sequences[0].seq)
    
    # 逐位生成consensus
    consensus = []
    for pos in range(seq_length):
        bases = [str(seq.seq[pos]) for seq in sequences]
        
        # 统计每个碱基的频率
        base_counts = Counter(bases)
        
        # 选择最常见的碱基
        most_common_base = base_counts.most_common(1)[0][0]
        consensus.append(most_common_base)
    
    return ''.join(consensus)

def optimize_large_cluster(aligned_file, output_file, max_size=500, sample_method='diverse', threads=4):
    """
    优化大cluster的consensus生成
    
    Args:
        aligned_file: 比对后的FASTA文件
        output_file: 输出consensus文件
        max_size: 采样后的最大大小
        sample_method: 采样方法
        threads: 线程数
    """
    logging.info(f"处理文件: {aligned_file}")
    
    try:
        sequences = list(SeqIO.parse(aligned_file, "fasta"))
        num_sequences = len(sequences)
        
        logging.info(f"序列数: {num_sequences}")
        
        # 根据序列数选择策略
        if num_sequences <= 500:
            logging.info("序列数<=500, 不需要优化,直接使用原文件")
            # 直接复制
            import shutil
            shutil.copy2(aligned_file, output_file)
            return
        
        elif num_sequences <= 1000:
            logging.info(f"序列数在500-1000之间, 使用{sample_method}采样到{max_size}条")
            # 中等采样
            sampled_sequences = representative_sampling(sequences, max_size, sample_method)
            SeqIO.write(sampled_sequences, output_file, "fasta")
            logging.info(f"采样完成: {len(sampled_sequences)} 条序列")
        
        else:
            logging.info(f"序列数>1000, 使用分层聚类+采样策略")
            
            # 先尝试分层聚类
            success = hierarchical_clustering_consensus(aligned_file, output_file, max_size, threads)
            
            if not success:
                # 如果分层聚类失败,使用激进采样
                logging.info(f"回退到{sample_method}采样策略")
                sampled_sequences = representative_sampling(sequences, max_size, sample_method)
                SeqIO.write(sampled_sequences, output_file, "fasta")
                logging.info(f"采样完成: {len(sampled_sequences)} 条序列")
    
    except Exception as e:
        logging.error(f"优化处理失败: {str(e)}")
        # 如果所有方法都失败,至少尝试随机采样
        try:
            sequences = list(SeqIO.parse(aligned_file, "fasta"))
            sampled_sequences = random.sample(sequences, min(max_size, len(sequences)))
            SeqIO.write(sampled_sequences, output_file, "fasta")
            logging.warning(f"使用随机采样作为最后手段: {len(sampled_sequences)} 条序列")
        except:
            logging.error("所有方法均失败,复制原文件")
            import shutil
            shutil.copy2(aligned_file, output_file)

def batch_optimize_clusters(input_dir, output_dir, max_size=500, sample_method='diverse', threads=4):
    """
    批量优化目录中的所有cluster文件
    """
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith('.fasta')]
    logging.info(f"找到 {len(fasta_files)} 个文件")
    
    for fasta_file in sorted(fasta_files):
        input_path = os.path.join(input_dir, fasta_file)
        output_path = os.path.join(output_dir, fasta_file)
        
        logging.info(f"\n处理: {fasta_file}")
        optimize_large_cluster(input_path, output_path, max_size, sample_method, threads)

def main():
    parser = argparse.ArgumentParser(description="优化大cluster的consensus生成")
    parser.add_argument("--input_file", help="输入aligned文件")
    parser.add_argument("--output_file", help="输出优化后的文件")
    parser.add_argument("--input_dir", help="输入目录(批量处理)")
    parser.add_argument("--output_dir", help="输出目录(批量处理)")
    parser.add_argument("--max_size", type=int, default=500,
                       help="采样后的最大序列数 (默认500)")
    parser.add_argument("--sample_method", choices=['random', 'diverse', 'stratified'],
                       default='diverse',
                       help="采样方法 (默认diverse)")
    parser.add_argument("--threads", type=int, default=4, help="线程数")
    
    args = parser.parse_args()
    
    if args.input_file:
        if not args.output_file:
            parser.error("单文件处理需要指定 --output_file")
        
        if not os.path.exists(args.input_file):
            logging.error(f"输入文件不存在: {args.input_file}")
            sys.exit(1)
        
        optimize_large_cluster(
            args.input_file,
            args.output_file,
            args.max_size,
            args.sample_method,
            args.threads
        )
    
    elif args.input_dir:
        if not args.output_dir:
            parser.error("批量处理需要指定 --output_dir")
        
        if not os.path.isdir(args.input_dir):
            logging.error(f"输入目录不存在: {args.input_dir}")
            sys.exit(1)
        
        batch_optimize_clusters(
            args.input_dir,
            args.output_dir,
            args.max_size,
            args.sample_method,
            args.threads
        )
    
    else:
        parser.error("必须指定 --input_file 或 --input_dir")

if __name__ == "__main__":
    main()
