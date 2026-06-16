#!/usr/bin/env python
# -*- coding: utf-8 -*-
# time: 2024/4/19 10:26
# author: Huyangzhirong
# file: amplicon_depth_stat.py
# -*- coding: utf-8 -*-

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
    with open(f"{bed}","r") as BED:
        dic_bed = {}
        bed_info = {}
        for line in BED:
            lines = line.strip("\n").split("\t")
            dic_bed[f"{lines[3]}_reads_depth"] = 0
            dic_bed[f"{lines[3]}_sc_depth"] = 0
            dic_bed.setdefault(f"{lines[3]}_insert_size","")
            bed_info[lines[3]] = lines
    return dic_bed, bed_info

def stat_reads(args):
    # 统计reads的函数
    # args: 包含df, chrom, pos, reads_le, strand的元组
    df, chrom, pos, reads_le, strand = args
    if strand == "+":
        intervals = df[(df['chrom'] == chrom) & ((df['start']-10) <= pos) & (df['end'] >= (pos+reads_le))]
    else:
        intervals = df[(df['chrom'] == chrom)  & ((df['end']+10) >= (pos+reads_le)) & (df['start'] <= pos)]
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
    sc_length=sum(length for op, length in read.cigartuples if op == 4)
    strand = "-" if read.is_reverse else "+"
    if read.is_paired and read.is_read1:
        insert_size = abs(read.template_length)
    else:
        insert_size = None
    return stat_reads((df, chrom, position, read_length, strand)),sc_length,insert_size

def stat_bam_split(i, dic_bed,df,bam,num_workers):
    # 统计每个bam文件分割部分中的reads
    split_reads = split_bam_in_memory(bam, num_workers)
    for read in split_reads[i]:
        result=process_read(read, df)
        if result:
            bed_result = result[0]
            sc_length = result[1]
            insert_size=result[2]
            if bed_result and bed_result != "NA":
                if insert_size:
                    insert_size_str=str(insert_size)+"|"
                    dic_bed[f"{bed_result}_insert_size"]+=insert_size_str
                dic_bed[f"{bed_result}_reads_depth"]+=1
                if sc_length>=80:
                    dic_bed[f"{bed_result}_sc_depth"]+=1

    return

def get_insert_size_gradient(bed_length,insert_size_list):
    insert_size_list.sort()
    insert_size_list_len=len(insert_size_list)
    gradient_list=[]
    for i in range(5):
        gradient_list.append(len([x for x in insert_size_list if x>=i*bed_length/5 and x<=(i+1)*bed_length/5])/insert_size_list_len)
    return gradient_list


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
        head=['chrom', "start", "end", "ID","depth","length","softclip_depth","softclip_depth_radio","insert_size_0%~20%","insert_size_20%~40%","insert_size_40%~60%","insert_size_60%~80%","insert_size_80%~100%\n"]
        OUT.write("\t".join(head))
        for key in bed_info:
            bed_length=int(bed_info[key][2])-int(bed_info[key][1])
            insert_size_list=shared_dict[f"{key}_insert_size"]
            if insert_size_list!="":
                insert_size_list=list(map(int,insert_size_list.split("|")[:-1]))
                insert_size_gradient=get_insert_size_gradient(bed_length,insert_size_list)
                insert_size_gradient_str="\t".join([f"{x:.2f}" for x in insert_size_gradient])
            else:
                insert_size_gradient_str="NA\tNA\tNA\tNA\tNA"
            softclip_depth_radio = "{:.2f}".format(shared_dict[f"{key}_sc_depth"] / shared_dict[f"{key}_reads_depth"]) if shared_dict[f"{key}_reads_depth"] != 0 else "NA"
            text = "\t".join(bed_info[key]) + "\t" + str(shared_dict[f"{key}_reads_depth"]) +"\t"+str(bed_length)+"\t"+str(shared_dict[f"{key}_sc_depth"]) +"\t"+str(softclip_depth_radio)+"\t"+insert_size_gradient_str+"\n"
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
