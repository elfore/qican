import os
import pandas as pd

def create_depth_matrix(project_name, input_file, output_dir):
    """
    将简洁版的表格转置为矩阵格式：
    行：区域 (Chrom, Start, End, Target)
    列：各样本的 SE1 和 SE2 深度值
    """
    if not os.path.exists(input_file):
        print(f"[{project_name}] Error: Input file {input_file} not found.")
        return

    print(f"\n=== Transposing Project: {project_name} ===")
    df = pd.read_csv(input_file, sep='\t')

    # 1. 提取 Depth 矩阵
    # 我们需要将 Depth_SE1 和 Depth_SE2 展开
    # 首先处理 SE1
    df_se1 = df[['SampleID', 'Chrom', 'Start', 'End', 'Target', 'Depth_SE1']].copy()
    df_se1['SampleID'] = df_se1['SampleID'] + "_SE1"
    df_se1.rename(columns={'Depth_SE1': 'Depth'}, inplace=True)

    # 然后处理 SE2
    df_se2 = df[['SampleID', 'Chrom', 'Start', 'End', 'Target', 'Depth_SE2']].copy()
    df_se2['SampleID'] = df_se2['SampleID'] + "_SE2"
    df_se2.rename(columns={'Depth_SE2': 'Depth'}, inplace=True)

    # 合并后透视
    df_combined = pd.concat([df_se1, df_se2])
    
    # 使用 pivot 将 SampleID 变为列
    matrix = df_combined.pivot(
        index=['Chrom', 'Start', 'End', 'Target'], 
        columns='SampleID', 
        values='Depth'
    ).reset_index()

    # 排序：按染色体和位置
    matrix.sort_values(by=['Chrom', 'Start'], inplace=True)

    # 保存矩阵
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"SKII15275_{project_name}_depth_matrix.txt")
    matrix.to_csv(out_path, sep='\t', index=False)
    
    print(f"[{project_name}] Matrix created: {out_path}")
    print(f"[{project_name}] Matrix shape: {matrix.shape[0]} rows x {matrix.shape[1]} columns")

# 配置
base_path = "/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/PE"

for proj in ['PA', 'PB', 'PD']:
    input_f = f"{base_path}/{proj}/SKII15275_{proj}_concise_depth.txt"
    create_depth_matrix(proj, input_f, f"{base_path}/{proj}")
