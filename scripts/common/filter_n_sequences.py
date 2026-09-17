#!/usr/bin/env python3
"""
Filter sequences with high N content and gaps
This script filters out sequences where the proportion of 'N' bases or '-' gaps exceeds 10%
"""

import argparse
import sys
from Bio import SeqIO
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_n_and_gap_percentage(sequence):
    """
    计算序列中N碱基和gap的百分比
    
    Args:
        sequence (str): 序列字符串
        
    Returns:
        tuple: (N碱基百分比, gap百分比, 总百分比)
    """
    if not sequence:
        return 0.0, 0.0, 0.0
    
    # 转换为大写并计算N和gap的数量
    seq_upper = sequence.upper()
    n_count = seq_upper.count('N')
    gap_count = seq_upper.count('-')
    total_length = len(sequence)
    
    if total_length == 0:
        return 0.0, 0.0, 0.0
    
    n_percentage = (n_count / total_length) * 100
    gap_percentage = (gap_count / total_length) * 100
    total_percentage = n_percentage + gap_percentage
    
    return n_percentage, gap_percentage, total_percentage

def filter_sequences(input_file, output_file, threshold=10.0):
    """
    过滤序列，移除N碱基和gap比例超过阈值的序列
    
    Args:
        input_file (str): 输入FASTA文件路径
        output_file (str): 输出FASTA文件路径
        threshold (float): N碱基和gap总比例阈值（百分比）
        
    Returns:
        tuple: (保留的序列数量, 过滤掉的序列数量)
    """
    try:
        # 读取输入文件
        sequences = list(SeqIO.parse(input_file, "fasta"))
        logging.info(f"读取了 {len(sequences)} 条序列")
        
        # 过滤序列
        filtered_sequences = []
        removed_count = 0
        
        for record in sequences:
            n_percentage, gap_percentage, total_percentage = calculate_n_and_gap_percentage(str(record.seq))
            
            if total_percentage <= threshold:
                filtered_sequences.append(record)
                logging.debug(f"保留序列 {record.id}: N比例 = {n_percentage:.2f}%, gap比例 = {gap_percentage:.2f}%, 总计 = {total_percentage:.2f}%")
            else:
                removed_count += 1
                logging.info(f"过滤序列 {record.id}: N比例 = {n_percentage:.2f}%, gap比例 = {gap_percentage:.2f}%, 总计 = {total_percentage:.2f}% (超过阈值 {threshold}%)")
        
        # 保存过滤后的序列
        SeqIO.write(filtered_sequences, output_file, "fasta")
        
        logging.info(f"过滤完成: 保留 {len(filtered_sequences)} 条序列，过滤掉 {removed_count} 条序列")
        logging.info(f"过滤后的序列已保存到: {output_file}")
        
        return len(filtered_sequences), removed_count
        
    except Exception as e:
        logging.error(f"过滤序列时发生错误: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Filter sequences with high N content and gaps")
    parser.add_argument("--input_file", required=True, help="Input FASTA file")
    parser.add_argument("--output_file", required=True, help="Output FASTA file")
    parser.add_argument("--threshold", type=float, default=10.0, 
                       help="N content and gap threshold percentage (default: 10.0)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # 验证输入文件
    if not args.input_file:
        logging.error("输入文件路径不能为空")
        sys.exit(1)
    
    # 执行过滤
    try:
        kept_count, removed_count = filter_sequences(
            args.input_file, 
            args.output_file, 
            args.threshold
        )
        
        # 输出统计信息
        print(f"过滤统计:")
        print(f"  输入序列数: {kept_count + removed_count}")
        print(f"  保留序列数: {kept_count}")
        print(f"  过滤序列数: {removed_count}")
        print(f"  N碱基和gap阈值: {args.threshold}%")
        
    except Exception as e:
        logging.error(f"程序执行失败: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 