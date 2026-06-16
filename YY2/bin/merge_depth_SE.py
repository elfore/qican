import os
import pandas as pd
import glob
import argparse

# usage:/mnt/gpfs1/Users/yangjinxurong/software/miniconda3/envs/nextflow/bin/python3 /mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/bin/merge_depth_SE.py -skid SKII15275 
# 适用于SE数据，用SE流程执行的结果，结果目录中包含result/depth和result/depth_MAPQ, 
# 最终report中使用的是depth_MAPQ的结果

def merge_project_concise(project_name, skid, se1_base, se2_base, output_dir):
    """
    将 SE_1 和 SE_2 数据合并到一行，使表格更简洁
    """
    def load_data(base_path, suffix):
        if not os.path.exists(base_path):
            return pd.DataFrame()
        
        d_dir = os.path.join(base_path, "result/depth") if os.path.exists(os.path.join(base_path, "result/depth")) else base_path
        m_dir = os.path.join(base_path, "result/depth_MAPQ") if os.path.exists(os.path.join(base_path, "result/depth_MAPQ")) else base_path
        
        d_files = glob.glob(os.path.join(d_dir, "*.reads_depth.txt"))
        all_df_list = []
        
        for df_path in d_files:
            sample_name = os.path.basename(df_path).replace(".reads_depth.txt", "")
            mf_path = os.path.join(m_dir, f"{sample_name}.reads_depth_MAPQ.txt")
            if not os.path.exists(mf_path):
                mf_path = os.path.join(m_dir, f"{sample_name}.reads_depth.txt")

            try:
                # 读取 Depth
                with open(df_path, 'r') as f:
                    first = f.readline()
                if first.startswith("chrom") or "start" in first:
                    df = pd.read_csv(df_path, sep='\t')
                    df.rename(columns={'chrom':'Chrom','start':'Start','end':'End','ID':'Target','depth':'Depth'}, inplace=True)
                else:
                    df = pd.read_csv(df_path, sep='\t', header=None).iloc[:, :5]
                    df.columns = ['Chrom', 'Start', 'End', 'Target', 'Depth']
                
                # 读取 MAPQ
                if os.path.exists(mf_path):
                    m_df = pd.read_csv(mf_path, sep='\t', header=None if not first.startswith("chrom") else 0)
                    df['Depth_MAPQ'] = m_df['depth'] if 'depth' in m_df.columns else m_df.iloc[:, 4]
                else:
                    df['Depth_MAPQ'] = '.'
                
                df['SampleID'] = sample_name
                all_df_list.append(df[['SampleID', 'Chrom', 'Start', 'End', 'Target', 'Depth', 'Depth_MAPQ']])
            except Exception as e:
                print(f"[{project_name}] Error loading {sample_name}: {e}")
        
        return pd.concat(all_df_list, ignore_index=True) if all_df_list else pd.DataFrame()

    print(f"\n=== Processing Project: {project_name} ===")
    df1 = load_data(se1_base, "SE_1")
    df2 = load_data(se2_base, "SE_2")

    if df1.empty or df2.empty:
        print(f"[{project_name}] Missing data for one or both SE directories.")
        return

    # 合并两个 DataFrame
    # 按 SampleID 和 坐标/Target 合并
    merged = pd.merge(
        df1, df2, 
        on=['SampleID', 'Chrom', 'Start', 'End', 'Target'], 
        suffixes=('_SE1', '_SE2')
    )

    # 调整列顺序
    cols = ['SampleID', 'Chrom', 'Start', 'End', 'Target', 'Depth_SE1', 'Depth_SE2', 'Depth_MAPQ_SE1', 'Depth_MAPQ_SE2']
    merged = merged[cols]

    # 排序
    merged.sort_values(by=['SampleID', 'Chrom', 'Start'], inplace=True)

    # 保存
    os.makedirs(output_dir, exist_ok=True)
    out_all = os.path.join(output_dir, f"{skid}_{project_name}_depth.txt")
    merged.to_csv(out_all, sep='\t', index=False)

    # 低深度过滤 (SE1 或 SE2 任一深度 < 20)
    merged['Depth_SE1'] = pd.to_numeric(merged['Depth_SE1'], errors='coerce')
    merged['Depth_SE2'] = pd.to_numeric(merged['Depth_SE2'], errors='coerce')
    lt20 = merged[(merged['Depth_SE1'] < 20) | (merged['Depth_SE2'] < 20)]
    
    out_lt20 = os.path.join(output_dir, f"{skid}_{project_name}_lt20.txt")
    lt20.to_csv(out_lt20, sep='\t', index=False)

    # 提取未过滤mappind rate的列生成汇报表格
    result_df = transpose_data(merged)
    output_file = os.path.join(output_dir, f"{skid}_{project_name}_depth_report.txt")
    result_df.to_csv(output_file, sep='\t', index=False)

    print(f"[{project_name}] Finished. Output: {output_file}")

def transpose_data(df):
    df_subset = df[['SampleID', 'Chrom', 'Start', 'End', 'Target', 'Depth_MAPQ_SE1', 'Depth_MAPQ_SE2']].copy()
    df_subset.columns = ['SampleID', 'Chrom', 'Start', 'End', 'Target', 'Depth_SE1', 'Depth_SE2']

    # 透视表：将 SampleID 展开为列
    pivot_df = df_subset.pivot_table(
        index=['Chrom', 'Start', 'End', 'Target'],
        columns='SampleID',
        values=['Depth_SE1', 'Depth_SE2'],
        aggfunc='first'
    )

    # 调整多级索引顺序，先按 SampleID 排序，再按 SE1/SE2 排序
    # 这样可以确保同一个样品的 SE1 和 SE2 紧挨在一起
    pivot_df = pivot_df.swaplevel(0, 1, axis=1).sort_index(axis=1)

    # 重新整理列名
    pivot_df.columns = [f"{sample}_{var}" for sample, var in pivot_df.columns]

    # 重置索引
    result_df = pivot_df.reset_index()

    # 生成示例（前5行，展示前10列，现在应该能看到同一个样品的 SE1 和 SE2 了）
    print("\n--- 结果示例 (前10列) ---")
    print(result_df.iloc[:5, :10].to_string(index=False))

    return result_df
    # 保存完整结果

# 配置
base_in = "/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150"
base_out = "/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/SE"

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-skid", required=True, help="skid")
    parse.add_argument("-single", help="SE", default=True)
    args = parse.parse_args()
    skid = args.skid
    for proj in ['PA', 'PB', 'PD']:
        merge_project_concise(proj, skid, f"{base_in}/{skid}_SE_1/{proj}", f"{base_in}/{skid}_SE_2/{proj}", f"{base_out}/{proj}")

if __name__ == "__main__":
    main()

