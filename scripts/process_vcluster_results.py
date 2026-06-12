#!/usr/bin/env python3
"""
Process vcluster clustering results and generate cluster information.

This script parses the vcluster output files and generates:
1. A cluster information CSV file with sequence-to-cluster mappings
2. A FASTA file with sequences organized by clusters
"""

import sys
import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import os

def parse_vcluster_results(clusters_file, ids_file):
    """
    Parse vcluster output files to extract clustering information.
    
    Args:
        clusters_file (str): Path to the clusters.tsv file
        ids_file (str): Path to the ani.ids.tsv file
        
    Returns:
        dict: Dictionary mapping sequence IDs to cluster IDs
    """
    cluster_mapping = {}
    
    # 读取序列ID映射文件 (vcluster格式: id\tseq_len\tno_parts)
    ids_df = pd.read_csv(ids_file, sep='\t')
    
    # 检查文件格式并相应处理
    if 'idx' in ids_df.columns and 'object' in ids_df.columns:
        # 旧格式: idx\tobject
        id_to_name = dict(zip(ids_df['idx'], ids_df['object']))
    elif 'id' in ids_df.columns:
        # 新格式: id\tseq_len\tno_parts
        # 创建索引到序列名的映射
        id_to_name = {i: row['id'] for i, row in ids_df.iterrows()}
    else:
        raise ValueError(f"Unsupported ids file format. Expected columns: 'idx'/'object' or 'id'. Found: {ids_df.columns.tolist()}")
    
    # 读取聚类结果文件
    clusters_df = pd.read_csv(clusters_file, sep='\t')
    
    # 检查聚类文件格式
    if 'object' in clusters_df.columns and 'cluster' in clusters_df.columns:
        # 新格式: object\tcluster
        for _, row in clusters_df.iterrows():
            sequence_name = row['object']
            cluster_id = row['cluster']
            cluster_mapping[sequence_name] = f"cluster_{cluster_id}"
    elif 'cluster' in clusters_df.columns and 'members' in clusters_df.columns:
        # 旧格式: cluster\tmembers
        for _, row in clusters_df.iterrows():
            cluster_id = row['cluster']
            # 获取该聚类中的所有序列ID
            if 'members' in row:
                member_indices = [int(x) for x in str(row['members']).split(',') if x.strip()]
                for idx in member_indices:
                    if idx in id_to_name:
                        sequence_name = id_to_name[idx]
                        cluster_mapping[sequence_name] = f"cluster_{cluster_id}"
    else:
        raise ValueError(f"Unsupported clusters file format. Expected columns: 'object'/'cluster' or 'cluster'/'members'. Found: {clusters_df.columns.tolist()}")
    
    return cluster_mapping

def create_cluster_info(cluster_mapping, aligned_fasta):
    """
    Create cluster information DataFrame.
    
    Args:
        cluster_mapping (dict): Sequence ID to cluster ID mapping
        aligned_fasta (str): Path to aligned FASTA file
        
    Returns:
        pd.DataFrame: Cluster information with sequence and cluster details
    """
    cluster_data = []
    
    # 读取序列以获取额外信息
    sequences = {}
    for record in SeqIO.parse(aligned_fasta, "fasta"):
        sequences[record.id] = str(record.seq)
    
    # 创建聚类信息
    for seq_id, cluster_id in cluster_mapping.items():
        if seq_id in sequences:
            seq_length = len(sequences[seq_id].replace('-', ''))  # 不含gap的长度
            cluster_data.append({
                'sequence_id': seq_id,
                'cluster_id': cluster_id,
                'sequence_length': seq_length,
                'aligned_sequence': sequences[seq_id]
            })
    
    return pd.DataFrame(cluster_data)

def write_clustered_sequences(cluster_mapping, aligned_fasta, output_file):
    """
    Write sequences organized by clusters to FASTA file.
    
    Args:
        cluster_mapping (dict): Sequence ID to cluster ID mapping
        aligned_fasta (str): Path to aligned FASTA file
        output_file (str): Output FASTA file path
    """
    # 按聚类分组序列
    clusters = defaultdict(list)
    for record in SeqIO.parse(aligned_fasta, "fasta"):
        if record.id in cluster_mapping:
            cluster_id = cluster_mapping[record.id]
            clusters[cluster_id].append(record)
    
    # 写入聚类序列，保持原始序列名称
    with open(output_file, 'w') as f:
        for cluster_id, records in clusters.items():
            for record in records:
                # 保持原始序列名称，不添加聚类前缀
                f.write(f">{record.id}\n{record.seq}\n")

def extract_genotype_from_filename(filename):
    """
    Extract genotype information from filename.
    
    Args:
        filename (str): Input filename
        
    Returns:
        str: Extracted genotype or "Unknown"
    """
    import re
    # 尝试从文件名中提取基因型信息
    genotype_match = re.search(r'(GII\.\d+)', filename)
    if genotype_match:
        return genotype_match.group(1)
    return "Unknown"

def process_vcluster_results(clusters_file, ids_file, aligned_fasta, cluster_output, sequences_output):
    """
    Main function to process vcluster clustering results.
    
    Args:
        clusters_file (str): Path to the clusters.tsv file
        ids_file (str): Path to the ani.ids.tsv file
        aligned_fasta (str): Path to aligned FASTA file
        cluster_output (str): Output CSV file for cluster information
        sequences_output (str): Output FASTA file for clustered sequences
    """
    print(f"Parsing vcluster clustering results from {clusters_file}...")
    
    # 解析vcluster结果
    cluster_mapping = parse_vcluster_results(clusters_file, ids_file)
    print(f"Found {len(cluster_mapping)} sequences in {len(set(cluster_mapping.values()))} clusters")
    
    # 创建聚类信息
    cluster_info = create_cluster_info(cluster_mapping, aligned_fasta)
    
    # 保存聚类信息
    cluster_info.to_csv(cluster_output, index=False)
    print(f"Cluster information saved to {cluster_output}")
    
    # 写入聚类序列
    write_clustered_sequences(cluster_mapping, aligned_fasta, sequences_output)
    print(f"Clustered sequences saved to {sequences_output}")
    
    # 打印统计摘要
    cluster_counts = cluster_info['cluster_id'].value_counts()
    print(f"\nCluster summary:")
    print(f"Total clusters: {len(cluster_counts)}")
    print(f"Largest cluster: {cluster_counts.max()} sequences")
    print(f"Smallest cluster: {cluster_counts.min()} sequences")
    print(f"Average cluster size: {cluster_counts.mean():.2f} sequences")
    
    # 分析基因型分布
    cluster_info['genotype'] = cluster_info['sequence_id'].str.extract(r'(GII\.\d+)')
    genotype_dist = cluster_info.groupby(['cluster_id', 'genotype']).size().unstack(fill_value=0)
    print(f"\nGenotype distribution across clusters:")
    print(genotype_dist.head(10))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process vcluster clustering results.")
    parser.add_argument("--clusters_file", type=str, required=True, help="Path to the clusters.tsv file from vcluster")
    parser.add_argument("--ids_file", type=str, required=True, help="Path to the ani.ids.tsv file from vcluster")
    parser.add_argument("--aligned_fasta", type=str, required=True, help="Path to the aligned FASTA file")
    parser.add_argument("--cluster_output", type=str, required=True, help="Output CSV file for cluster information")
    parser.add_argument("--sequences_output", type=str, required=True, help="Output FASTA file for clustered sequences")
    
    args = parser.parse_args()
    
    process_vcluster_results(args.clusters_file, args.ids_file, args.aligned_fasta, 
                           args.cluster_output, args.sequences_output) 