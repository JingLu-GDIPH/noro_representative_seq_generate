#!/usr/bin/env python3
"""
Split clustered sequences into individual FASTA files per cluster.

This script takes the cluster information and clustered sequences,
then creates separate FASTA files for each cluster with genotype information in filenames.
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def extract_genotype(sequence_id):
    """Extract genotype from sequence_id (e.g., GII.4_MZ546180 -> GII.4)"""
    return sequence_id.split('_')[0]

def get_cluster_genotype(cluster_sequences):
    """Get genotype information for a cluster"""
    genotypes = set()
    for seq_id in cluster_sequences:
        genotype = extract_genotype(seq_id)
        genotypes.add(genotype)
    
    # Sort genotypes for consistent naming
    sorted_genotypes = sorted(genotypes)
    if len(sorted_genotypes) == 1:
        return sorted_genotypes[0]
    else:
        # 使用简化的命名方式，避免文件名过长
        if len(sorted_genotypes) <= 2:
            return '_'.join(sorted_genotypes)
        else:
            # 如果基因型太多，使用混合命名
            return f"mixed_{len(sorted_genotypes)}_genotypes"

def create_safe_filename(cluster_id, genotype, max_length=100):
    """Create a safe filename with length limit"""
    base_name = f"cluster__{cluster_id}__{genotype}"
    if len(base_name) <= max_length:
        return f"{base_name}.fasta"
    else:
        # 如果太长，截断基因型部分
        available_length = max_length - len(f"cluster__{cluster_id}__") - len(".fasta")
        if available_length > 0:
            truncated_genotype = genotype[:available_length]
            return f"cluster__{cluster_id}__{truncated_genotype}.fasta"
        else:
            # 如果还是太长，使用最简单的命名
            return f"cluster__{cluster_id}.fasta"

def split_clusters(cluster_info_file, sequences_file, outdir):
    """
    Split clustered sequences into individual FASTA files per cluster.
    
    Args:
        cluster_info_file (str): Path to cluster info CSV file
        sequences_file (str): Path to clustered sequences FASTA file
        outdir (str): Output directory for cluster FASTA files
    """
    # Create output directory
    os.makedirs(outdir, exist_ok=True)
    
    # Read cluster information
    print(f"Reading cluster information from {cluster_info_file}...")
    df = pd.read_csv(cluster_info_file)
    
    # Group sequences by cluster
    cluster_groups = df.groupby('cluster_id')['sequence_id'].apply(list).to_dict()
    print(f"Found {len(cluster_groups)} clusters")
    
    # Read all sequences
    print(f"Reading sequences from {sequences_file}...")
    seq_dict = {}
    for record in SeqIO.parse(sequences_file, "fasta"):
        seq_dict[record.id] = record
    
    # Create FASTA file for each cluster
    cluster_files = []
    for cluster_id, seq_ids in cluster_groups.items():
        # Get genotype information for this cluster. Newer genotype-pair clustering
        # writes this explicitly; older vclust-based tables fall back to ID parsing.
        if "genotype_pair" in df.columns:
            genotype_values = sorted(
                value for value in df.loc[df["cluster_id"] == cluster_id, "genotype_pair"].dropna().unique()
            )
            genotype = genotype_values[0] if len(genotype_values) == 1 else "_".join(genotype_values)
        else:
            genotype = get_cluster_genotype(seq_ids)
        
        # Create filename with genotype information
        filename = create_safe_filename(cluster_id, genotype)
        out_fasta = os.path.join(outdir, filename)
        records = [seq_dict[sid] for sid in seq_ids if sid in seq_dict]
        
        if records:
            SeqIO.write(records, out_fasta, "fasta")
            cluster_files.append(out_fasta)
            print(f"Cluster {cluster_id} ({genotype}): {len(records)} sequences -> {out_fasta}")
        else:
            print(f"Warning: No sequences found for cluster {cluster_id}")
    
    print(f"\nSuccessfully created {len(cluster_files)} cluster FASTA files in {outdir}")
    return cluster_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split clustered FASTA into per-cluster files with genotype information.")
    parser.add_argument("--cluster_info", required=True, help="CSV with sequence_id and cluster_id columns")
    parser.add_argument("--sequences", required=True, help="Clustered sequences FASTA")
    parser.add_argument("--outdir", required=True, help="Output directory for per-cluster FASTA files")
    
    args = parser.parse_args()
    
    split_clusters(args.cluster_info, args.sequences, args.outdir)
