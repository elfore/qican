import os
import gzip
import pandas as pd
import glob
import argparse

# usage:/mnt/gpfs1/Users/yangjinxurong/software/miniconda3/envs/nextflow/bin/python3 /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/bin/merge_depth_PE.py -skid SKII15275

EXPECTED_COLS = ['SampleID', 'Chrom', 'Start', 'End', 'Target', 'Depth']


def merge_project_concise(project_name, skid, pe_base, output_dir):
    """
    Merge depth files for one project and output both the full depth report and
    a report-shaped table containing rows where at least one sample has depth <20.
    """
    def find_depth_files(base_path):
        if not os.path.exists(base_path):
            return []

        candidate_dirs = [
            os.path.join(base_path, "result/depth"),
            os.path.join(base_path, "result2/depth"),
            base_path,
        ]
        for d_dir in candidate_dirs:
            if not os.path.exists(d_dir):
                continue
            d_files = sorted(glob.glob(os.path.join(d_dir, "*.reads_depth.txt")))
            if d_files:
                return d_files

        # Some PB runs did not publish reads_depth files into result/depth.
        work_files = sorted(glob.glob(os.path.join(base_path, "work/*/*/*.reads_depth.txt")))
        if work_files:
            return work_files

        # Last-resort PB hotspot depth output. The first five columns still map
        # to Chrom, Start, End, Target, Depth.
        for d_dir in [
            os.path.join(base_path, "result/hotspot_depth"),
            os.path.join(base_path, "result2/hotspot_depth"),
        ]:
            if not os.path.exists(d_dir):
                continue
            d_files = sorted(glob.glob(os.path.join(d_dir, "*.regions.bed.gz")))
            if d_files:
                return d_files

        return []

    def open_text(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

    def sample_name_from_path(path):
        name = os.path.basename(path)
        for suffix in [".reads_depth.txt", ".regions.bed.gz"]:
            if name.endswith(suffix):
                return name[:-len(suffix)]
        return os.path.splitext(name)[0]

    def normalize_columns(df, path):
        rename_map = {
            '#chrom': 'Chrom', 'chrom': 'Chrom', 'chr': 'Chrom',
            'start': 'Start', 'end': 'End',
            'ID': 'Target', 'id': 'Target', 'name': 'Target', 'target': 'Target',
            'depth': 'Depth', 'Depth': 'Depth', 'mean': 'Depth', 'mean_depth': 'Depth',
        }
        df.rename(columns=rename_map, inplace=True)
        missing = [col for col in ['Chrom', 'Start', 'End', 'Target', 'Depth'] if col not in df.columns]
        if missing:
            raise ValueError(f"missing columns {missing} in {path}")
        return df[['Chrom', 'Start', 'End', 'Target', 'Depth']]

    def read_depth_file(df_path):
        sample_name = sample_name_from_path(df_path)
        with open_text(df_path) as f:
            first = f.readline()

        first_lower = first.lower()
        has_header = first_lower.startswith(("chrom", "#chrom")) or "start" in first_lower
        if has_header:
            df = pd.read_csv(df_path, sep='\t', compression='infer')
            df = normalize_columns(df, df_path)
        else:
            df = pd.read_csv(df_path, sep='\t', header=None, compression='infer')
            if df.shape[1] < 5:
                raise ValueError(f"expected at least 5 columns, got {df.shape[1]} in {df_path}")
            df = df.iloc[:, :5]
            df.columns = ['Chrom', 'Start', 'End', 'Target', 'Depth']

        df['SampleID'] = sample_name
        return df[EXPECTED_COLS]

    def load_data(base_path):
        d_files = find_depth_files(base_path)
        if not d_files:
            print(f"[{project_name}] Warning: no depth files found under {base_path}")
            return pd.DataFrame(columns=EXPECTED_COLS)

        all_df_list = []
        print(f"[{project_name}] Loading {len(d_files)} depth files")
        for df_path in d_files:
            sample_name = sample_name_from_path(df_path)
            try:
                all_df_list.append(read_depth_file(df_path))
            except Exception as e:
                print(f"[{project_name}] Error loading {sample_name}: {e}")

        return pd.concat(all_df_list, ignore_index=True) if all_df_list else pd.DataFrame(columns=EXPECTED_COLS)

    print(f"\n=== Processing Project: {project_name} ===")

    merged = load_data(pe_base)
    merged = merged[EXPECTED_COLS]

    os.makedirs(output_dir, exist_ok=True)
    out_all = os.path.join(output_dir, f"{skid}_{project_name}_depth.txt")
    output_file = os.path.join(output_dir, f"{skid}_{project_name}_depth_report.txt")
    out_lt20 = os.path.join(output_dir, f"{skid}_{project_name}_lt20_report.txt")

    if merged.empty:
        merged.to_csv(out_all, sep='\t', index=False)
        pd.DataFrame(columns=['Chrom', 'Start', 'End', 'Target']).to_csv(output_file, sep='\t', index=False)
        pd.DataFrame(columns=['Chrom', 'Start', 'End', 'Target']).to_csv(out_lt20, sep='\t', index=False)
        print(f"[{project_name}] No valid depth data. Empty outputs written to {output_dir}")
        return

    merged['Depth'] = pd.to_numeric(merged['Depth'], errors='coerce')
    merged.sort_values(by=['SampleID', 'Chrom', 'Start'], inplace=True)
    merged.to_csv(out_all, sep='\t', index=False)

    result_df = transpose_data(merged)
    result_df.to_csv(output_file, sep='\t', index=False)

    depth_cols = [col for col in result_df.columns if col.endswith('_Depth')]
    lt20_mask = result_df[depth_cols].apply(pd.to_numeric, errors='coerce').lt(20).any(axis=1)
    result_lt20 = result_df.loc[lt20_mask].copy()
    result_lt20.to_csv(out_lt20, sep='\t', index=False)

    print(f"[{project_name}] Finished. Output: {output_file}")
    print(f"[{project_name}] Low-depth report rows: {len(result_lt20)}. Output: {out_lt20}")


def transpose_data(df):
    df_subset = df[EXPECTED_COLS].copy()

    # 透视表：将 SampleID 展开为列
    pivot_df = df_subset.pivot_table(
        index=['Chrom', 'Start', 'End', 'Target'],
        columns='SampleID',
        values=['Depth'],
        aggfunc='first'
    )

    # 调整多级索引顺序，按 SampleID 排序
    pivot_df = pivot_df.swaplevel(0, 1, axis=1).sort_index(axis=1)

    # 重新整理列名
    pivot_df.columns = [f"{sample}_{var}" for sample, var in pivot_df.columns]

    # 重置索引
    result_df = pivot_df.reset_index()

    # 生成示例（前5行，展示前10列）
    print("\n--- 结果示例 (前10列) ---")
    print(result_df.iloc[:5, :10].to_string(index=False))

    return result_df


# 配置
base_in = "/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150"
base_out = "/mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/PE"


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-skid", required=True, help="skid")
    parse.add_argument("-single", help="SE", default=True)
    args = parse.parse_args()
    skid = args.skid
    for proj in ['PA', 'PB', 'PD']:
        merge_project_concise(proj, skid, f"{base_in}/{skid}/{proj}", f"{base_out}/{proj}")


if __name__ == "__main__":
    main()
