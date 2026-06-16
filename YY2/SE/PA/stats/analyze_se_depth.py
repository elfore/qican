import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

def analyze():
    input_file = "../SKII15275_PA_depth_processed.txt"
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    df = pd.read_csv(input_file, sep='\t')
    se1_cols = [c for c in df.columns if c.endswith('_Depth_SE1')]
    se2_cols = [c for c in df.columns if c.endswith('_Depth_SE2')]
    samples = [c.replace('_Depth_SE1', '') for c in se1_cols]

    # 1. 样本维度
    sample_stats = []
    for s in samples:
        se1 = df[f"{s}_Depth_SE1"]
        se2 = df[f"{s}_Depth_SE2"]
        corr, _ = stats.pearsonr(se1, se2) if se1.sum() > 0 and se2.sum() > 0 else (np.nan, np.nan)
        sample_stats.append({
            'Sample': s,
            'Mean_SE1': se1.mean(), 'Mean_SE2': se2.mean(),
            'Pearson_R': corr,
            'Ratio': se1.mean() / se2.mean() if se2.mean() != 0 else np.nan
        })
    pd.DataFrame(sample_stats).to_csv("summary_by_sample.csv", index=False)

    # 2. Target 维度
    target_info = df[['Chrom', 'Start', 'End', 'Target']].copy()
    target_info['Avg_SE1'] = df[se1_cols].mean(axis=1)
    target_info['Avg_SE2'] = df[se2_cols].mean(axis=1)
    target_info['Ratio'] = target_info['Avg_SE1'] / target_info['Avg_SE2']
    
    val_se1 = df[se1_cols].values
    val_se2 = df[se2_cols].values
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = val_se1 / val_se2
    target_info['Ratio_Std'] = np.nanstd(ratios, axis=1)
    target_info.to_csv("summary_by_target.csv", index=False)

    print("Python 统计汇总完成。")

if __name__ == "__main__":
    analyze()
