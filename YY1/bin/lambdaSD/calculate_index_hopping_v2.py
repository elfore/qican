import sys
import pysam
import glob
import os
import pandas as pd

mapfile = sys.argv[1]
mapdict = {}
with open(mapfile) as f:
    for line in f:
        line = line.strip()
        linelist = line.split('\t')
        mapdict[linelist[1]] = linelist[0]

def count_read_align(bam, mapdict):
    bam_name = os.path.basename(bam).split('.')[0].replace('Read1-', '')
    bf = pysam.AlignmentFile(bam, 'rb')
    chr_length = dict(zip(bf.references,bf.lengths))
    count_correct = 0 
    count_not_mapped = 0 
    count_misassigned = {}  # 用于存储比对到错误的参考名的计数
    total_count = 0

    for r in bf.fetch(until_eof=True):
        if r.is_unmapped:
            count_not_mapped += 1 
        else:
            if r.is_duplicate or r.is_supplementary or r.is_secondary or r.is_qcfail:
                continue
            else:
                if r.is_forward:
                    if r.reference_start <=10:
                        total_count += 1
                        if r.reference_name == mapdict[bam_name]:
                            count_correct += 1
                        elif  r.reference_name != mapdict[bam_name]:
                            misassigned_ref = r.reference_name
                            for sample, ref in mapdict.items():
                                if misassigned_ref == ref:
                                    if sample not in count_misassigned:
                                        count_misassigned[sample] = 0
                                    count_misassigned[sample] += 1
                    else:
                        continue
                elif r.is_reverse:
                    if r.reference_end >= int(chr_length[r.reference_name])-10:
                        total_count += 1
                        if r.reference_name == mapdict[bam_name]:
                            count_correct += 1
                        elif  r.reference_name != mapdict[bam_name]:
                            misassigned_ref = r.reference_name
                            for sample, ref in mapdict.items():
                                if misassigned_ref == ref:
                                    if sample not in count_misassigned:
                                        count_misassigned[sample] = 0
                                    count_misassigned[sample] += 1
                    else:
                        continue

    return total_count, count_correct, count_not_mapped, count_misassigned

bam_dir = sys.argv[2]
bams = glob.glob(bam_dir + '/*.bam')
results = []
misassigned_counts = {}  # 用于存储所有样本的错误比对计数

# 初始化 misassigned_counts
for bam in bams:
    bam_name = os.path.basename(bam).split('.')[0].replace('Read1-', '')
    total_count, count_correct, count_not_mapped, count_misassigned = count_read_align(bam, mapdict)
    
    # 保存每个样本的计数
    results.append({
        'Sample': bam_name,
        'Total_Reads': total_count,
        'Correctly_Mapped': count_correct,
        'Misassigned_Counts': sum(count_misassigned.values()),
        'Unmapped': count_not_mapped
    })

    # 更新 misassigned_counts
    misassigned_counts[bam_name] = count_misassigned

# 创建结果 DataFrame
df = pd.DataFrame(results)
df = df.sort_values(by='Sample')
# 创建样本之间的矩阵
sample_names = sorted(list(df['Sample']))
matrix = pd.DataFrame(0, index=sample_names, columns=sample_names)

# 填充矩阵
for sample, misassigned in misassigned_counts.items():
    for misassigned_sample, count in misassigned.items():
        if misassigned_sample in sample_names:
            matrix.loc[sample, misassigned_sample] = count

# 将对角线填充为正确比对的计数
for index, row in df.iterrows():
    matrix.loc[row['Sample'], row['Sample']] = row['Correctly_Mapped']

# 保存结果
df.to_csv('read_align_count.xls', index=False, sep="\t")
matrix.to_csv('misassigned_counts_matrix.xls', sep="\t")
