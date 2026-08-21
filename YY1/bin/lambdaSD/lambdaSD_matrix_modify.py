# -*- coding: utf-8 -*-
# time: 2024/10/23 13:41
# file: lambdaSD_matrix_modify.py
import sys

import pandas as pd

df = pd.read_csv(sys.argv[1], sep="\t", index_col=0)

# 获取矩阵的形状
rows, cols = df.shape

# 创建一个新的DataFrame用于存储结果
results = pd.DataFrame(columns=['sampleID','Total_Reads', 'Correctly_Mapped ', 'Misassigned_Counts'])
# 遍历每一行
for i in range(len(df)):
    # 检查对角线周围的3个单元格
    for j in range(-3, 4):
        if j != 0:  # 不包括对角线本身
            if 0 <= i + j < len(df):
                if df.iloc[i + j, i] != 0:
                    df.iloc[i + j, i] = 0  # 改写为0


Correctly_Mapped = df.values.diagonal()
Total_Reads = df.sum(axis=1)
Misassigned_Counts = Total_Reads - Correctly_Mapped

results = pd.DataFrame({
    'sampleID': df.index,
    'Total_Reads':Total_Reads,
    'Correctly_Mapped':Correctly_Mapped,
    'Misassigned_Counts':Misassigned_Counts
})

# 输出结果到新的CSV文件
df.to_csv("misassigned_counts_matrix_filter.xls", sep="\t")
results.to_csv('read_align_count_filter.xls', sep="\t", index=False)

# 汇总输出过滤版 index hopping rate
total_reads_sum = int(results['Total_Reads'].sum())
correctly_mapped_sum = int(results['Correctly_Mapped'].sum())
misassigned_sum = int(results['Misassigned_Counts'].sum())
index_hopping_rate = misassigned_sum / total_reads_sum if total_reads_sum else 0

rate_df = pd.DataFrame([{
    'Total_Reads': total_reads_sum,
    'Correctly_Mapped': correctly_mapped_sum,
    'Misassigned_Counts': misassigned_sum,
    'Index_Hopping_Rate': index_hopping_rate,
    'Index_Hopping_Rate_Percent': index_hopping_rate * 100,
    'Rate_Type': 'filtered_diagonal_neighbor_pm3'
}])
rate_df.to_csv('index_hopping_rate.xls', sep="\t", index=False)
