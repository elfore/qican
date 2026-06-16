import argparse
import pysam
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import os
import time
import psutil
from joblib import Parallel, delayed
import multiprocessing

def get_bed_dic(bed):
    # 读取bed文件，返回两个字典
    # dic_bed: 存储bed文件中的第四列作为键，值初始化为0
    # bed_info: 存储bed文件中的每一行信息，以第四列作为键
    with open(f"{bed}", "r") as BED:
        dic_bed = {}
        bed_info = {}
        for line in BED:
            lines = line.strip("\n").split("\t")
            dic_bed[lines[3]] = 0
            bed_info[lines[3]] = lines
    return dic_bed, bed_info

def stat_reads(args):
    # 统计reads的函数
    # args: 包含df, chrom, pos, reads_le, strand的元组
    df, chrom, pos, reads_le, strand = args
    if strand == "+":
        intervals = df[(df['chrom'] == chrom) & ((df['start']+15) >= pos) & ((df['start']-15) <= pos) & ((df['end'] + 20)>= (pos+reads_le))]
    else:
        intervals = df[(df['chrom'] == chrom) & ((df['end']+15) >= (pos+reads_le)) & ((df['end']-15) <= (pos+reads_le)) & ((df['start'] - 20)<= pos)]
    intervals_list = intervals.values.tolist()
    if len(intervals_list) == 1:
        return intervals_list[0][3]
    else:
        return "NA"

def split_bam_in_memory(bam_path, num_splits):
    # 将bam文件分割为多个部分，存储在内存中
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    reads = [read for read in bamfile.fetch()]
    bamfile.close()
    
    split_size = len(reads) // num_splits

    split_reads = [reads[i:i + split_size] for i in range(0, len(reads), split_size)]

    if len(split_reads) > num_splits:
        last_reads = split_reads.pop()
        split_reads[-1].extend(last_reads)
    
    return split_reads

def process_read(read, df):
    # 处理每个read，返回统计结果
    if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary or read.is_qcfail:
        return None
    chrom = read.reference_name
    position = read.reference_start
    read_length = sum(length for op, length in read.cigartuples if op == 0 or op == 2)
    strand = "-" if read.is_reverse else "+"
    return stat_reads((df, chrom, position, read_length, strand))

def stat_bam_split(i, dic_bed,df,bam,num_workers):
    # 统计每个bam文件分割部分中的reads
    split_reads = split_bam_in_memory(bam, num_workers)
    for read in split_reads[i]:
        result = process_read(read, df)
        if result and result != "NA":
            dic_bed[result] += 1
    return 


def get_bam(bed, bam, outfile, num_workers):
    # 主函数，获取bed和bam文件，进行统计
    start_time = time.time()
    process = psutil.Process(os.getpid())
    start_mem = process.memory_info().rss
    dic_bed, bed_info = get_bed_dic(bed)
    manager = multiprocessing.Manager()
    shared_dict = manager.dict()
    shared_dict.update(dic_bed)
    df = pd.read_csv(bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'primer'], dtype={'chrom': str, 'start': int, 'end': int, 'primer': str})
    # 使用并行处理统计每个bam文件分割部分中的reads
    Parallel(n_jobs=num_workers)(delayed(stat_bam_split)(i, shared_dict,df,bam,num_workers) for i in range(num_workers))
    with open(f"{outfile}", "w") as OUT:
        for key in dic_bed:
            text = "\t".join(bed_info[key]) + "\t" + str(shared_dict[key]) + "\n"
            OUT.write(text)
    end_time = time.time()
    end_mem = process.memory_info().rss
    print(f"Execution time: {end_time - start_time:.2f} seconds")
    print(f"Memory used: {(end_mem - start_mem) / (1024**2):.2f} MB")


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-bed", required=True, help="bed文件路径")
    parse.add_argument("-bam", required=True, help="bam文件路径")
    parse.add_argument("-outfile", required=True, help="输出文件路径")
    parse.add_argument("-workers", type=int, default=8, help="并行处理的工作进程数")
    args = parse.parse_args()
    get_bam(args.bed, args.bam, args.outfile, args.workers)

if __name__ == "__main__":
    main()
