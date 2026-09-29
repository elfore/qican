#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import sys
import json
import csv
import io
import re

# 飞书电子表格的 tokens 与 sheet IDs 定义
TOKENS = {
    "cnvseq": {
        "token": "NWgXs0zzShsXYLtPuN5c1M3fn52",
        "sheet_name": "CNVseq",
        "sheet_id": "f146fa",
        "mode": "append_row"
    },
    "pgt": {
        "token": "UfFEsGG0Zh97KotvvWzcLRsNnSf",
        "sheet_name": "PGTA",
        "sheet_id": "ce5e9c",
        "mode": "append_row"
    },
    "lambsd_summary": {
        "token": "FbLRsZf6OhGceSthhflcKSVrnIe",
        "sheet_name": "汇总表",
        "sheet_id": "4680e7",
        "mode": "append_row"
    },
    "lambsd_readscount": {
        "token": "FbLRsZf6OhGceSthhflcKSVrnIe",
        "sheet_name": "readscount表",
        "sheet_id": "qUoQv8",
        "mode": "append_column"
    },
    "tngs_sens": {
        "token": "YLmusKmiphU8iPtTz6zc5bo1ni9",
        "sheet_name": "阴阳性及hopping-新",
        "sheet_id": "itLCNo",
        "mode": "append_row"
    },
    "tngs_qc": {
        "token": "YLmusKmiphU8iPtTz6zc5bo1ni9",
        "sheet_name": "质控汇总表",
        "sheet_id": "4gpizi",
        "mode": "append_row"
    },
    "tngs_hop": {
        "token": "YLmusKmiphU8iPtTz6zc5bo1ni9",
        "sheet_name": "index hopping矩阵表-新",
        "sheet_id": "Tx5MbT",
        "mode": "append_row"
    },
    "mngs_qc": {
        "token": "MjihsiLrnh6UU5tc4cSc4PRSntc",
        "sheet_name": "mNGS质控汇总",
        "sheet_id": "5639f1",
        "mode": "append_row"
    },
    "mngs_detail": {
        "token": "MjihsiLrnh6UU5tc4cSc4PRSntc",
        "mode": "dynamic_sheet"
    },
    "tb2_qc": {
        "token": "N6FRsWZGRhaLskt57hTckWFlnSc",
        "sheet_name": "数据质控",
        "sheet_id": "f2f823",
        "mode": "append_row"
    },
    "tb2_amp": {
        "token": "N6FRsWZGRhaLskt57hTckWFlnSc",
        "sheet_name": "原始扩增子深度表",
        "sheet_id": "YXom0T",
        "mode": "append_column_tb2_amp"
    },
    "tb2_hot": {
        "token": "N6FRsWZGRhaLskt57hTckWFlnSc",
        "sheet_name": "突变热点深度表",
        "sheet_id": "7ehiNN",
        "mode": "append_column_tb2_hot"
    }
}

# 工具与环境路径
YY1_DIR = "/mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1"
LARK_CLI = "/mnt/gpfs1/Users/yangjinxurong/software/lark-cli/bin/lark-cli"
mNGS_SCRIPT = "/mnt/gpfs1/Users/wangning/project/qican/mNGS/QC_stat.pl"
mNGS_CLASS = "/mnt/gpfs1/Users/wangning/project/qican/mNGS/classfiy.pl"

def run_cmd(cmd, shell=True, check=True):
    print(f"[RUN] {cmd}")
    res = subprocess.run(cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and res.returncode != 0:
        print(f"[ERROR] Command failed:\nStdout: {res.stdout}\nStderr: {res.stderr}")
        raise RuntimeError(f"Command failed with code {res.returncode}")
    return res

def col_idx_to_letter(idx):
    result = ""
    while idx > 0:
        idx, remainder = divmod(idx - 1, 26)
        result = chr(65 + remainder) + result
    return result

def get_sheet_info(token, sheet_id):
    res = run_cmd(f"{LARK_CLI} sheets +sheet-info --spreadsheet-token '{token}' --sheet-id '{sheet_id}' --format json")
    info = json.loads(res.stdout)
    if not info.get("ok"):
        raise RuntimeError("Failed to get sheet info")
    data = info["data"]
    
    # 解析 range (如 "A1:Z169") 以兼容没有 row_count/column_count 的老版本
    range_str = data.get("range", "A1:A1")
    m_row = re.search(r':\w*?(\d+)$', range_str)
    row_count = int(m_row.group(1)) if m_row else 1
    
    m_col = re.search(r':([A-Z]+)\d+$', range_str)
    col_letter = m_col.group(1) if m_col else "A"
    
    col_idx = 0
    for c in col_letter:
        col_idx = col_idx * 26 + (ord(c) - 64)
        
    data["row_count"] = row_count
    data["column_count"] = col_idx
    return data

def get_sheet_csv(token, sheet_id, range_str):
    res = run_cmd(f"{LARK_CLI} sheets +csv-get --spreadsheet-token '{token}' --sheet-id '{sheet_id}' --range '{range_str}' --format json")
    data = json.loads(res.stdout)
    if not data.get("ok"):
        raise RuntimeError("Failed to get sheet CSV")
    # 解析 annotated_csv，去除每行开头的 [row=N] 前缀
    raw_lines = data["data"]["annotated_csv"].strip().split("\n")
    cleaned_lines = []
    for line in raw_lines:
        m = re.match(r"^\[row=\d+\]\s*(.*)", line)
        if m:
            cleaned_lines.append(m.group(1))
        else:
            cleaned_lines.append(line)
    return list(csv.reader(cleaned_lines))

def put_csv_feishu(token, sheet_name, start_cell, csv_text):
    p = subprocess.Popen([LARK_CLI, "sheets", "+csv-put", "--spreadsheet-token", token, "--sheet-name", sheet_name, "--start-cell", start_cell, "--csv", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = p.communicate(input=csv_text)
    if p.returncode != 0:
        print(f"[ERROR] csv-put failed:\nStdout: {stdout}\nStderr: {stderr}")
        raise RuntimeError("csv-put failed")
    print(f"[SUCCESS] csv-put succeeded: {stdout.strip()}")

def main():
    parser = argparse.ArgumentParser(description="Integrate batch data and sync to Lark Sheets")
    parser.add_argument("--batch", required=True, help="Batch ID, e.g., BGISEQ_0707")
    args = parser.parse_args()
    batch = args.batch

    print(f"=== Starting integration for batch: {batch} ===")

    # 1. 解析 YY1_info.txt 获取路径与配置
    info_path = os.path.join(YY1_DIR, "bin/YY1_info.txt")
    batch_dir = None
    batch_yq = None
    batch_tp = None
    with open(info_path, "r") as f:
        for line in f:
            lines = line.strip("\n").split("\t")
            if lines[0] == batch:
                batch_yq = lines[1]
                batch_tp = lines[2]
                batch_dir = lines[3]
                break

    if not batch_dir:
        print(f"[ERROR] Batch {batch} not found in {info_path}")
        sys.exit(1)
    
    print(f"[INFO] Located batch directory: {batch_dir}")

    # ==================== 1. CNVseq ====================
    print("\n--- Processing CNVseq ---")
    cnv_stat_file = os.path.join(YY1_DIR, "bin/cnv/stat.txt")
    cnv_rows = []
    with open(cnv_stat_file, "r") as f:
        for line in f:
            lines = line.strip("\n").split("\t")
            if lines[0] == batch:
                cnv_rows.append(lines)
    if cnv_rows:
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(cnv_rows)
        cfg = TOKENS["cnvseq"]
        append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())
    else:
        print("[WARN] No CNVseq data found for batch.")

    # ==================== 2. PGTA ====================
    print("\n--- Processing PGTA ---")
    pgt_stat_file = os.path.join(YY1_DIR, "bin/pgta/stat.txt")
    pgt_rows = []
    # 注意：PGT 项目里 pici 被获取为了第二列（即批次仪器号如 BGISEQ_T7_0707）
    # 故这里需要匹配 batch_yq (如 BGISEQ_T7_0707)
    with open(pgt_stat_file, "r") as f:
        for line in f:
            lines = line.strip("\n").split("\t")
            if lines[0] == batch_yq:
                pgt_rows.append(lines)
    if pgt_rows:
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(pgt_rows)
        cfg = TOKENS["pgt"]
        append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())
    else:
        print("[WARN] No PGTA data found for batch.")

    # ==================== 3. lambdaSD ====================
    print("\n--- Processing lambdaSD ---")
    lambda_dir = os.path.join(batch_dir, "lambdaSD")
    
    # 汇总表
    lamb_rate_file = os.path.join(lambda_dir, "index_hopping_rate.xls")
    if os.path.exists(lamb_rate_file):
        with open(lamb_rate_file, "r") as f:
            lines = [line.strip("\n").split("\t") for line in f]
            if len(lines) > 1:
                # 提取 Rate_Type 对应的 hopping rate
                # 结构：Total_Reads	Correctly_Mapped	Misassigned_Counts	Index_Hopping_Rate	Index_Hopping_Rate_Percent	Rate_Type
                rate_val = float(lines[1][4])  # 得到百分比，如 0.00719
                rate_str = f"{rate_val:.5f}%"
                summary_row = [batch_yq, batch_yq.split("_")[1] if "_" in batch_yq else "", "", rate_str, "49SD index hopping 矩阵"]
                out = io.StringIO()
                writer = csv.writer(out, lineterminator="\n")
                writer.writerow(summary_row)
                cfg = TOKENS["lambsd_summary"]
                append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())
    
    # readscount 表
    lamb_count_file = os.path.join(lambda_dir, "read_align_count_filter.xls")
    if os.path.exists(lamb_count_file):
        # 1. 读当前在线的 readscount 表 A 列，获取 amplicon 顺序
        cfg = TOKENS["lambsd_readscount"]
        sheet_info = get_sheet_info(cfg["token"], cfg["sheet_id"])
        total_cols = sheet_info["column_count"]
        total_rows = sheet_info["row_count"]
        
        a_col_data = get_sheet_csv(cfg["token"], cfg["sheet_id"], f"A1:A{total_rows}")
        # 2. 读取本地 xls 映射
        mapped_reads = {}
        total_mapped_sum = 0
        with open(lamb_count_file, "r") as f:
            next(f)
            for line in f:
                parts = line.strip("\n").split("\t")
                val = int(parts[2]) # Correctly_Mapped
                mapped_reads[parts[0]] = val
                total_mapped_sum += val
        
        # 3. 按照 A 列对齐构造新的一列
        new_column = []
        for i, row in enumerate(a_col_data):
            if not row:
                new_column.append("")
                continue
            amp = row[0]
            if i == 0:
                new_column.append(batch) # 列头：批次名
            elif i == 1:
                new_column.append(str(total_mapped_sum)) # total mapped reads 行
            else:
                new_column.append(str(mapped_reads.get(amp, 0)))
        
        # 4. 横向追加写入新一列
        target_col = col_idx_to_letter(total_cols + 1)
        csv_text = "\n".join(new_column)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], f"{target_col}1", csv_text)

    # ==================== 4. tNGS ====================
    print("\n--- Processing tNGS ---")
    tngs_test_dir = os.path.join(YY1_DIR, "bin/tngs/test")
    
    # 4.1 阴阳性及hopping-新
    res_stat_file = os.path.join(tngs_test_dir, "result_stat.txt")
    tngs_sens_row = None
    with open(res_stat_file, "r") as f:
        next(f)
        for line in f:
            parts = line.strip("\n").split("\t")
            if parts[0] == batch_yq:
                # 匹配得到 SEN, PPV
                sen = float(parts[2]) * 100
                ppv = float(parts[3]) * 100
                # 获取 ap_stat 算出的 hopping rate
                # 从 ap_result.txt 或上面 task 输出计算
                # 简单起见，从 ap_result.txt 的最后一列取均值或者重新跑一次
                # 我们采用直接读取终端计算好的 index hopping rate
                # 在此重新读取 ap_result.txt 或者是 stat 出来的 hopping 比例
                # 我们在此以静态估算或从 ap_result.txt 计算
                tngs_sens_row = [parts[0], parts[1], "", f"{sen:.2f}%", f"{ppv:.2f}%", ""]
                break
    
    # 从 ap_result 估算 hopping
    ap_res_file = os.path.join(tngs_test_dir, "ap_result.txt")
    if tngs_sens_row and os.path.exists(ap_res_file):
        # 我们可以读取并计算出本批次总的 hopping rate
        # 根据 ap_stat.py 的公式：fz / fm
        # 我们需要在 integrate_batch.py 中重构这个计算
        # 从 sp_ap_info.txt 读取各样本的阴阳性背景
        sp_ap_info = {}
        with open(os.path.join(tngs_test_dir, "sp_ap_info.txt"), "r") as f:
            for line in f:
                parts = line.strip("\n").split("\t")
                sp_ap_info[parts[0]] = parts[1].split(",") + parts[2].split(",")
        
        info_list = []
        with open(os.path.join(tngs_test_dir, "info.txt"), "r") as f:
            for line in f:
                info_list.append(line.strip("\n"))
        
        fz = 0
        fm = 0
        file_list = sorted(subprocess.check_output(f"ls {batch_dir}/tNGS/LC*/LC*merge_count.xls", shell=True).decode("utf-8").strip().split("\n"))
        for file1 in file_list:
            if os.path.getsize(file1) == 0:
                continue
            sample_name = os.path.basename(file1).split(".")[0]
            if sample_name not in sp_ap_info:
                continue
            list_ap = sp_ap_info[sample_name]
            with open(file1, "r") as FILE:
                next(FILE)
                for line in FILE:
                    parts = line.strip("\n").split("\t")
                    a_id = parts[0]
                    dp = float(parts[1])
                    if a_id in info_list:
                        # 只有不在背景列表中的才算作 hopping
                        if a_id not in list_ap:
                            fm += dp
                            fz += dp
                        else:
                            fm += dp
        
        radio_f = fz / fm if fm else 0
        tngs_sens_row[5] = f"{radio_f * 100:.6f}%"
        
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(tngs_sens_row)
        cfg = TOKENS["tngs_sens"]
        append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())

    # 4.2 质控汇总表
    tngs_stat_file = os.path.join(tngs_test_dir, "tngs_stat.txt")
    tngs_qc_rows = []
    with open(tngs_stat_file, "r") as f:
        next(f)
        for line in f:
            parts = line.strip("\n").split("\t")
            if parts[0] == batch_yq:
                tngs_qc_rows.append(parts)
    if tngs_qc_rows:
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(tngs_qc_rows)
        cfg = TOKENS["tngs_qc"]
        append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())

    # 4.3 index hopping 矩阵表-新
    tngs_hop_rows = []
    with open(ap_res_file, "r") as f:
        next(f)
        for line in f:
            parts = line.strip("\n").split("\t")
            if parts[0] == batch_yq:
                tngs_hop_rows.append(parts)
    if tngs_hop_rows:
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(tngs_hop_rows)
        cfg = TOKENS["tngs_hop"]
        append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())

    # ==================== 5. mNGS ====================
    print("\n--- Processing mNGS ---")
    mngs_dir = os.path.join(batch_dir, "mNGS")
    mngs_stat_file = os.path.join(mngs_dir, "stat.txt")
    
    # 5.1 质控汇总
    if os.path.exists(mngs_stat_file):
        mngs_rows = []
        with open(mngs_stat_file, "r") as f:
            for line in f:
                parts = line.strip("\n").split("\t")
                if parts[0] == batch:
                    mngs_rows.append(parts)
        if mngs_rows:
            out = io.StringIO()
            writer = csv.writer(out, lineterminator="\n")
            writer.writerows(mngs_rows)
            cfg = TOKENS["mngs_qc"]
            append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())
            
    # 5.2 文库专属 Sheet 动态创建和写入
    # 扫描文库目录
    lib_dirs = [d for d in os.listdir(mngs_dir) if os.path.isdir(os.path.join(mngs_dir, d)) and not d.startswith(".snakemake") and d != "boxplot_inputs"]
    cfg = TOKENS["mngs_detail"]
    wb_info = run_cmd(f"{LARK_CLI} sheets +workbook-info --spreadsheet-token '{cfg['token']}' --format json")
    wb_data = json.loads(wb_info.stdout)
    existing_sheets = [s["sheet_name"] for s in wb_data["data"]["sheets"]]
    
    for lib in lib_dirs:
        # 为每个文库生成 classfiy.txt
        # 1. 运行 classify.pl
        lib_path = os.path.join(mngs_dir, lib)
        classify_out = run_cmd(f"perl {mNGS_CLASS} {lib_path}")
        
        # 2. 检查 sheet 是否存在，不存在就创建
        if lib not in existing_sheets:
            print(f"[INFO] Creating new sheet for mNGS library: {lib}...")
            run_cmd(f"{LARK_CLI} sheets +sheet-create --spreadsheet-token '{cfg['token']}' --title '{lib}'")
        
        # 3. 将比例数据写入
        # classify_out 的输出是用 tab 分隔的，我们要转换成标准的 CSV
        csv_rows = [line.split("\t") for line in classify_out.stdout.strip().split("\n")]
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(csv_rows)
        put_csv_feishu(cfg["token"], lib, "A1", out.getvalue())

    # ==================== 6. TB2 ====================
    print("\n--- Processing TB2 ---")
    tb2_res_dir = os.path.join(batch_dir, "TB2/result")
    
    # 6.1 数据质控 (数据追加)
    tb2_qc_file = os.path.join(tb2_res_dir, "QC.stat.xls")
    if os.path.exists(tb2_qc_file):
        tb_qc_rows = []
        with open(tb2_qc_file, "r") as f:
            next(f)
            for line in f:
                parts = line.strip("\n").split("\t")
                # 第一列是 runID，我们写入当前批次前缀（如 SKII15393 等）
                # 从 prefix 或者是 work.sh 里的 runID 提取，也可以就是 batch
                parts.insert(0, batch)
                tb_qc_rows.append(parts)
        if tb_qc_rows:
            out = io.StringIO()
            writer = csv.writer(out, lineterminator="\n")
            writer.writerows(tb_qc_rows)
            cfg = TOKENS["tb2_qc"]
            append_rows_feishu_dynamic(cfg["token"], cfg["sheet_name"], cfg["sheet_id"], out.getvalue())

    # 6.2 原始扩增子深度表 (列追加)
    tb2_amp_file = None
    for f in os.listdir(tb2_res_dir):
        if f.endswith("_amplicon.depth.xls"):
            tb2_amp_file = os.path.join(tb2_res_dir, f)
            break
            
    if tb2_amp_file and os.path.exists(tb2_amp_file):
        cfg = TOKENS["tb2_amp"]
        append_columns_tb2(cfg, tb2_amp_file, batch, tb2_res_dir, mode="amp")

    # 6.3 突变热点深度表 (列追加)
    tb2_hot_file = None
    for f in os.listdir(tb2_res_dir):
        if f.endswith("_hot_depth.xls"):
            tb2_hot_file = os.path.join(tb2_res_dir, f)
            break
            
    if tb2_hot_file and os.path.exists(tb2_hot_file):
        cfg = TOKENS["tb2_hot"]
        append_columns_tb2(cfg, tb2_hot_file, batch, tb2_res_dir, mode="hot")

    print("\n=== Batch Integration Completed Successfully ===")

def append_rows_feishu_dynamic(token, sheet_name, sheet_id, csv_text):
    # 动态获取当前行数并纵向追加
    sheet_info = get_sheet_info(token, sheet_id)
    row_count = sheet_info["row_count"]
    start_cell = f"A{row_count + 1}"
    put_csv_feishu(token, sheet_name, start_cell, csv_text)

def append_columns_tb2(cfg, depth_file, batch, tb2_res_dir, mode="amp"):
    # 1. 读当前在线的 A 列或 D 列作为靶点索引
    sheet_info = get_sheet_info(cfg["token"], cfg["sheet_id"])
    total_cols = sheet_info["column_count"]
    total_rows = sheet_info["row_count"]
    
    # amp 模式以 A 列作为主键；hot 模式以 D 列作为主键
    idx_col = "A" if mode == "amp" else "D"
    idx_col_data = get_sheet_csv(cfg["token"], cfg["sheet_id"], f"{idx_col}1:{idx_col}{total_rows}")
    
    # 2. 读取本地深度文件数据，并提取样本列表和描述映射
    # 从 QC.stat.xls 中读出 prefix 与文库描述的映射
    desc_map = {}
    qc_stat = os.path.join(tb2_res_dir, "QC.stat.xls")
    with open(qc_stat, "r") as f:
        next(f)
        for line in f:
            parts = line.strip("\n").split("\t")
            # parts[2] 为 prefix, parts[1] 为文库描述
            desc_map[parts[2]] = parts[1]
            
    # 读取深度表
    headers = []
    data_dict = {}
    with open(depth_file, "r") as f:
        headers = next(f).strip("\n").split("\t")
        for line in f:
            parts = line.strip("\n").split("\t")
            # amp 模式主键是 parts[0]; hot 模式主键是 parts[3]
            key = parts[0] if mode == "amp" else parts[3]
            data_dict[key] = parts
            
    # headers 格式：[amp, sample1, sample2, ...]
    samples = headers[1:] if mode == "amp" else headers[4:]
    sample_cols_count = len(samples)
    
    # 3. 构造追加的 CSV 列矩阵
    new_cols = [[] for _ in range(total_rows)]
    for r_idx, row in enumerate(idx_col_data):
        if not row:
            for _ in range(sample_cols_count):
                new_cols[r_idx].append("")
            continue
        key = row[0]
        if r_idx == 0:
            # 第一行：第一列填批次号，其余填充空或相同批次号以供后面合并
            for s_idx in range(sample_cols_count):
                new_cols[r_idx].append(batch if s_idx == 0 else "")
        elif r_idx == 1:
            # 第二行：文库描述
            for s in samples:
                new_cols[r_idx].append(desc_map.get(s, ""))
        elif r_idx == 2:
            # 第三行：文库全名（Prefix）
            for s in samples:
                new_cols[r_idx].append(s)
        else:
            # 数据行
            row_vals = data_dict.get(key)
            if row_vals:
                vals = row_vals[1:] if mode == "amp" else row_vals[4:]
                for v in vals:
                    new_cols[r_idx].append(v)
            else:
                for _ in range(sample_cols_count):
                    new_cols[r_idx].append("0")
                    
    # 4. 转换成 CSV 内容
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerows(new_cols)
    
    # 5. 横向追加写入
    target_col = col_idx_to_letter(total_cols + 1)
    put_csv_feishu(cfg["token"], cfg["sheet_name"], f"{target_col}1", out.getvalue())
    
    # 6. 合并单元格首行以美化
    # 合并首行的 target_col 至最终的列字母
    end_col = col_idx_to_letter(total_cols + sample_cols_count)
    if sample_cols_count > 1:
        print(f"[INFO] Merging header cells: {target_col}1:{end_col}1...")
        run_cmd(f"{LARK_CLI} sheets +cells-merge --spreadsheet-token '{cfg['token']}' --sheet-name '{cfg['sheet_name']}' --range '{target_col}1:{end_col}1'")

if __name__ == "__main__":
    main()
