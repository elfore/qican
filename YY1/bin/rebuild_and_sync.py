#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import csv
import io
import re
import subprocess

# 飞书电子表格的 tokens 与 sheet IDs 定义
TOKENS = {
    "cnvseq": {
        "token": "NWgXs0zzShsXYLtPuN5c1M3fn52",
        "sheet_name": "CNVseq",
        "sheet_id": "f146fa",
        "header_rows": 1,
        "clean_range": "A2:Z1000",
        "start_cell": "A2"
    },
    "pgt": {
        "token": "UfFEsGG0Zh97KotvvWzcLRsNnSf",
        "sheet_name": "PGTA",
        "sheet_id": "ce5e9c",
        "header_rows": 2,
        "clean_range": "A3:Z1000",
        "start_cell": "A3"
    },
    "lambsd_summary": {
        "token": "FbLRsZf6OhGceSthhflcKSVrnIe",
        "sheet_name": "汇总表",
        "sheet_id": "4680e7",
        "header_rows": 1,
        "clean_range": "A2:Z1000",
        "start_cell": "A2"
    },
    "lambsd_readscount": {
        "token": "FbLRsZf6OhGceSthhflcKSVrnIe",
        "sheet_name": "readscount表",
        "sheet_id": "qUoQv8",
        "clean_range": "A1:AZ200",
        "start_cell": "A1"
    },
    "tngs_sens": {
        "token": "YLmusKmiphU8iPtTz6zc5bo1ni9",
        "sheet_name": "阴阳性及hopping-新",
        "sheet_id": "itLCNo",
        "header_rows": 1,
        "clean_range": "A2:Z1000",
        "start_cell": "A2"
    },
    "tngs_qc": {
        "token": "YLmusKmiphU8iPtTz6zc5bo1ni9",
        "sheet_name": "质控汇总表",
        "sheet_id": "4gpizi",
        "header_rows": 1,
        "clean_range": "A2:Z1000",
        "start_cell": "A2"
    },
    "tngs_hop": {
        "token": "YLmusKmiphU8iPtTz6zc5bo1ni9",
        "sheet_name": "index hopping矩阵表-新",
        "sheet_id": "Tx5MbT",
        "header_rows": 1,
        "clean_range": "A2:BF1000",
        "start_cell": "A2"
    },
    "mngs_qc": {
        "token": "MjihsiLrnh6UU5tc4cSc4PRSntc",
        "sheet_name": "mNGS质控汇总",
        "sheet_id": "5639f1",
        "header_rows": 2,
        "clean_range": "A3:Z1000",
        "start_cell": "A3"
    },
    "tb2_qc": {
        "token": "N6FRsWZGRhaLskt57hTckWFlnSc",
        "sheet_name": "数据质控",
        "sheet_id": "f2f823",
        "header_rows": 1,
        "clean_range": "A2:Z1000",
        "start_cell": "A2"
    },
    "tb2_amp": {
        "token": "N6FRsWZGRhaLskt57hTckWFlnSc",
        "sheet_name": "原始扩增子深度表",
        "sheet_id": "YXom0T",
        "clean_range": "A1:BZ300",
        "start_cell": "A1"
    },
    "tb2_hot": {
        "token": "N6FRsWZGRhaLskt57hTckWFlnSc",
        "sheet_name": "突变热点深度表",
        "sheet_id": "7ehiNN",
        "clean_range": "A1:BZ1000",
        "start_cell": "A1"
    }
}

YY1_DIR = "/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY1"
LARK_CLI = "/mnt/gpfs1/Users/yangjinxurong/software/lark-cli/bin/lark-cli"

def run_cmd(cmd, check=True):
    print(f"[RUN] {cmd}")
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and res.returncode != 0:
        print(f"[ERROR] Command failed with code {res.returncode}:\nStdout: {res.stdout}\nStderr: {res.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    return res

def clear_range(token, sheet_name, range_str):
    try:
        res = run_cmd(f"{LARK_CLI} sheets +sheet-info --spreadsheet-token '{token}' --sheet-name '{sheet_name}' --format json")
        info = json.loads(res.stdout)
        if info.get("ok"):
            for mc in info["data"].get("merged_cells", []):
                r_str = mc["range"]
                print(f"[INFO] Unmerging merged cell range before clear: {r_str}...")
                run_cmd(f"{LARK_CLI} sheets +cells-unmerge --spreadsheet-token '{token}' --sheet-name '{sheet_name}' --range '{r_str}'")
    except Exception as e:
        print(f"[WARN] Failed to unmerge cells for {sheet_name}: {e}")

    print(f"[INFO] Clearing online sheet {sheet_name} range {range_str}...")
    run_cmd(f"{LARK_CLI} sheets +cells-clear --spreadsheet-token '{token}' --sheet-name '{sheet_name}' --range '{range_str}' --yes", check=True)

def put_csv_feishu(token, sheet_name, start_cell, csv_text):
    p = subprocess.Popen([LARK_CLI, "sheets", "+csv-put", "--spreadsheet-token", token, "--sheet-name", sheet_name, "--start-cell", start_cell, "--csv", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = p.communicate(input=csv_text)
    if p.returncode != 0:
        print(f"[ERROR] csv-put failed:\nStdout: {stdout}\nStderr: {stderr}")
        raise RuntimeError("csv-put failed")
    print(f"[SUCCESS] csv-put succeeded: {stdout.strip()}")

def first_existing_path(*paths):
    """Prefer legacy layouts, with compatibility for current result layouts."""
    for path in paths:
        if os.path.exists(path):
            return path
    return paths[0]


def first_existing_dir(*paths):
    """Prefer legacy result directories, with a root-level TB2 fallback."""
    for path in paths:
        if os.path.isdir(path):
            return path
    return paths[0]


def col_idx_to_letter(idx):
    result = ""
    while idx > 0:
        idx, remainder = divmod(idx - 1, 26)
        result = chr(65 + remainder) + result
    return result

def get_sheet_csv(token, sheet_id, range_str):
    res = run_cmd(f"{LARK_CLI} sheets +csv-get --spreadsheet-token '{token}' --sheet-id '{sheet_id}' --range '{range_str}' --format json")
    data = json.loads(res.stdout)
    if not data.get("ok"):
        raise RuntimeError("Failed to get sheet CSV")
    raw_lines = data["data"]["annotated_csv"].strip().split("\n")
    cleaned_lines = []
    for line in raw_lines:
        m = re.match(r"^\[row=\d+\]\s*(.*)", line)
        if m:
            cleaned_lines.append(m.group(1))
        else:
            cleaned_lines.append(line)
    return list(csv.reader(cleaned_lines))

def main():
    print("=== Starting Rebuilding and Syncing All Sheets (Strictly Clean Format) ===")
    
    # 1. 读 YY1_info.txt 获取所有批次路径
    info_path = os.path.join(YY1_DIR, "bin/YY1_info.txt")
    batches = []
    with open(info_path, "r") as f:
        for line in f:
            if line.strip():
                parts = line.strip("\n").split("\t")
                batches.append({
                    "batch": parts[0],
                    "yq": parts[1],
                    "tp": parts[2],
                    "dir": parts[3]
                })

    # ==================== 1. CNVseq ====================
    print("\n--- Rebuilding CNVseq ---")
    cnv_file = os.path.join(YY1_DIR, "bin/cnv/stat.txt")
    if os.path.exists(cnv_file):
        rows = []
        with open(cnv_file, "r") as f:
            for line in f:
                if line.strip():
                    rows.append(line.strip("\n").split("\t"))
        if rows:
            cfg = TOKENS["cnvseq"]
            clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
            out = io.StringIO()
            writer = csv.writer(out, lineterminator="\n")
            writer.writerows(rows)
            put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # ==================== 2. PGTA ====================
    print("\n--- Rebuilding PGTA ---")
    pgt_file = os.path.join(YY1_DIR, "bin/pgta/stat.txt")
    if os.path.exists(pgt_file):
        rows = []
        with open(pgt_file, "r") as f:
            for line in f:
                if line.strip():
                    rows.append(line.strip("\n").split("\t"))
        if rows:
            cfg = TOKENS["pgt"]
            clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
            out = io.StringIO()
            writer = csv.writer(out, lineterminator="\n")
            writer.writerows(rows)
            put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # ==================== 3. lambdaSD ====================
    print("\n--- Rebuilding lambdaSD ---")
    # 3.1 汇总表
    summary_rows = []
    for b in batches:
        b_rate_file = first_existing_path(os.path.join(b["dir"], "lambdaSD/index_hopping_rate.xls"), os.path.join(b["dir"], "lambdaSD/results/index_hopping_rate.xls"))
        if os.path.exists(b_rate_file):
            with open(b_rate_file, "r") as f:
                lines = [line.strip("\n").split("\t") for line in f]
                if len(lines) > 1:
                    rate_val = float(lines[1][4])
                    rate_str = f"{rate_val:.5f}%"
                    summary_rows.append([b["yq"], b["yq"].split("_")[1] if "_" in b["yq"] else "", "", rate_str, "49SD index hopping 矩阵"])
    if summary_rows:
        cfg = TOKENS["lambsd_summary"]
        clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(summary_rows)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # 3.2 readscount表
    cfg = TOKENS["lambsd_readscount"]
    # 先在线获取 A 列的 amplicon 列表以对齐
    a_col_data = get_sheet_csv(cfg["token"], cfg["sheet_id"], "A1:A100")
    if a_col_data:
        # 清洗 A 列列表
        amplicon_list = [r[0] for r in a_col_data if r]
        
        # 扫描各个批次，重新构建整个二维矩阵
        cols_data = [[] for _ in range(len(amplicon_list))]
        # 第一列为 amplicon 名
        for i, amp in enumerate(amplicon_list):
            cols_data[i].append(amp)
            
        # 第二列为 ilmn定值，从在线原表的 B 列获取
        b_col_data = get_sheet_csv(cfg["token"], cfg["sheet_id"], f"B1:B{len(amplicon_list)}")
        for i, val in enumerate(b_col_data):
            cols_data[i].append(val[0] if val else "0")
            
        # 后续各列为各批次数据
        for b in batches:
            b_count_file = first_existing_path(os.path.join(b["dir"], "lambdaSD/read_align_count_filter.xls"), os.path.join(b["dir"], "lambdaSD/results/read_align_count_filter.xls"))
            if os.path.exists(b_count_file):
                mapped_reads = {}
                total_mapped_sum = 0
                with open(b_count_file, "r") as f:
                    next(f)
                    for line in f:
                        parts = line.strip("\n").split("\t")
                        val = int(parts[2])
                        mapped_reads[parts[0]] = val
                        total_mapped_sum += val
                
                # 填充该批次列
                for i, amp in enumerate(amplicon_list):
                    if i == 0:
                        cols_data[i].append(b["batch"]) # 批次名
                    elif i == 1:
                        cols_data[i].append(str(total_mapped_sum))
                    else:
                        cols_data[i].append(str(mapped_reads.get(amp, 0)))
                        
        cfg = TOKENS["lambsd_readscount"]
        clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(cols_data)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # ==================== 4. tNGS ====================
    print("\n--- Rebuilding tNGS ---")
    tngs_test_dir = os.path.join(YY1_DIR, "bin/tngs/test")
    
    # 4.1 阴阳性及hopping-新
    # 遍历所有批次，读取 result_stat.txt 并且结合估算 hopping 组装行
    tngs_sens_rows = []
    # 我们先解析 sp_ap_info.txt 与 info.txt
    sp_ap_info = {}
    with open(os.path.join(tngs_test_dir, "sp_ap_info.txt"), "r") as f:
        for line in f:
            parts = line.strip("\n").split("\t")
            sp_ap_info[parts[0]] = parts[1].split(",") + parts[2].split(",")
    info_list = []
    with open(os.path.join(tngs_test_dir, "info.txt"), "r") as f:
        for line in f:
            info_list.append(line.strip("\n"))

    # 从 result_stat.txt 获取批次敏感度
    sens_map = {}
    with open(os.path.join(tngs_test_dir, "result_stat.txt"), "r") as f:
        next(f)
        for line in f:
            parts = line.strip("\n").split("\t")
            sens_map[parts[0]] = parts
            
    for b in batches:
        if b["yq"] in sens_map:
            parts = sens_map[b["yq"]]
            sen = float(parts[2]) * 100
            ppv = float(parts[3]) * 100
            
            # 计算该批次的 hopping rate
            fz = 0
            fm = 0
            # 查找并计算该批次下的 hopping
            try:
                file_list = sorted(subprocess.check_output(f"ls {b['dir']}/tNGS/LC*/LC*merge_count.xls", shell=True).decode("utf-8").strip().split("\n"))
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
                            sub_parts = line.strip("\n").split("\t")
                            a_id = sub_parts[0]
                            dp = float(sub_parts[1])
                            if a_id in info_list:
                                if a_id not in list_ap:
                                    fm += dp
                                    fz += dp
                                else:
                                    fm += dp
                radio_f = fz / fm if fm else 0
                hop_str = f"{radio_f * 100:.6f}%"
            except Exception:
                hop_str = ""
                
            tngs_sens_rows.append([parts[0], parts[1], "", f"{sen:.2f}%", f"{ppv:.2f}%", hop_str])
            
    if tngs_sens_rows:
        cfg = TOKENS["tngs_sens"]
        clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(tngs_sens_rows)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # 4.2 质控汇总表
    tngs_qc_rows = []
    with open(os.path.join(tngs_test_dir, "tngs_stat.txt"), "r") as f:
        next(f)
        for line in f:
            if line.strip():
                tngs_qc_rows.append(line.strip("\n").split("\t"))
    if tngs_qc_rows:
        cfg = TOKENS["tngs_qc"]
        clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(tngs_qc_rows)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # 4.3 index hopping 矩阵表-新
    tngs_hop_rows = []
    with open(os.path.join(tngs_test_dir, "ap_result.txt"), "r") as f:
        next(f)
        for line in f:
            if line.strip():
                tngs_hop_rows.append(line.strip("\n").split("\t"))
    if tngs_hop_rows:
        cfg = TOKENS["tngs_hop"]
        clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(tngs_hop_rows)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # ==================== 5. mNGS ====================
    print("\n--- Rebuilding mNGS ---")
    mngs_rows = []
    for b in batches:
        b_mngs_stat = os.path.join(b["dir"], "mNGS/stat.txt")
        if os.path.exists(b_mngs_stat):
            with open(b_mngs_stat, "r") as f:
                for line in f:
                    if line.strip():
                        parts = line.strip("\n").split("\t")
                        # 确保第一列为正确的批次名
                        parts[0] = b["batch"]
                        mngs_rows.append(parts)
    if mngs_rows:
        cfg = TOKENS["mngs_qc"]
        clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(mngs_rows)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # ==================== 6. TB2 ====================
    print("\n--- Rebuilding TB2 ---")
    # 6.1 数据质控
    tb_qc_rows = []
    tb2_batches = [] # 记录含有 TB2 的批次
    for b in batches:
        b_tb_qc = first_existing_path(os.path.join(b["dir"], "TB2/result/QC.stat.xls"), os.path.join(b["dir"], "TB2/QC.stat.xls"))
        if os.path.exists(b_tb_qc):
            tb2_batches.append(b)
            with open(b_tb_qc, "r") as f:
                next(f)
                for line in f:
                    if line.strip():
                        parts = line.strip("\n").split("\t")
                        parts.insert(0, b["batch"])
                        tb_qc_rows.append(parts)
    if tb_qc_rows:
        cfg = TOKENS["tb2_qc"]
        clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerows(tb_qc_rows)
        put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())

    # 6.2 原始扩增子深度表 & 突变热点深度表 (全量重构矩阵)
    rebuild_tb2_depth_matrix(TOKENS["tb2_amp"], tb2_batches, mode="amp")
    rebuild_tb2_depth_matrix(TOKENS["tb2_hot"], tb2_batches, mode="hot")

    print("\n=== Rebuild and Sync Completed Successfully ===")

def rebuild_tb2_depth_matrix(cfg, tb2_batches, mode="amp"):
    # 1. 读当前在线的 A 列或 D 列作为靶点主键索引
    # 深度非常长，为了保险，读到 1000 行
    idx_col = "A" if mode == "amp" else "D"
    idx_col_data = get_sheet_csv(cfg["token"], cfg["sheet_id"], f"{idx_col}1:{idx_col}850")
    idx_list = [r[0] for r in idx_col_data if r]
    
    # 2. 依次读取每个包含 TB2 项目的批次的深度数据
    # 先组装所有的样本列
    # 每一列包含：第 1 行批次名，第 2 行文库描述，第 3 行 Prefix，第 4 行及以后是深度
    matrix = [[] for _ in range(len(idx_list))]
    
    # 填充主键及其他前置列
    if mode == "amp":
        # 只有第一列 amp
        for r_idx, key in enumerate(idx_list):
            matrix[r_idx].append(key)
    else:
        # hot 模式有四列：chr, start, end, mutation
        # 我们可以直接读取原表前四列并填入
        four_cols = get_sheet_csv(cfg["token"], cfg["sheet_id"], f"A1:D{len(idx_list)}")
        for r_idx, row in enumerate(four_cols):
            for c_idx in range(4):
                matrix[r_idx].append(row[c_idx] if c_idx < len(row) else "")

    # 合并单元格的范围列表，记录为 (start_col_letter, end_col_letter)
    merge_ranges = []
    current_col_idx = 2 if mode == "amp" else 5
    
    for b in tb2_batches:
        tb2_res_dir = first_existing_dir(os.path.join(b["dir"], "TB2/result"), os.path.join(b["dir"], "TB2"))
        depth_file = None
        suffix = "_amplicon.depth.xls" if mode == "amp" else "_hot_depth.xls"
        for f in os.listdir(tb2_res_dir):
            if f.endswith(suffix):
                depth_file = os.path.join(tb2_res_dir, f)
                break
                
        if not depth_file or not os.path.exists(depth_file):
            continue
            
        # 读描述
        desc_map = {}
        with open(os.path.join(tb2_res_dir, "QC.stat.xls"), "r") as f:
            next(f)
            for line in f:
                parts = line.strip("\n").split("\t")
                desc_map[parts[2]] = parts[1]
                
        # 读深度数据
        headers = []
        data_dict = {}
        with open(depth_file, "r") as f:
            headers = next(f).strip("\n").split("\t")
            for line in f:
                parts = line.strip("\n").split("\t")
                key = parts[0] if mode == "amp" else parts[3]
                data_dict[key] = parts
                
        samples = headers[1:] if mode == "amp" else headers[4:]
        sample_count = len(samples)
        if sample_count == 0:
            continue
            
        # 记录合并单元格的范围
        start_letter = col_idx_to_letter(current_col_idx)
        end_letter = col_idx_to_letter(current_col_idx + sample_count - 1)
        merge_ranges.append((start_letter, end_letter))
        current_col_idx += sample_count
        
        # 填充这一批次下的各个样本列
        for r_idx, key in enumerate(idx_list):
            if r_idx == 0:
                for s_idx in range(sample_count):
                    matrix[r_idx].append(b["batch"] if s_idx == 0 else "")
            elif r_idx == 1:
                # hot 模式的第二行在原飞书是“文库描述”吗？我们看前面 annotated_csv
                # 对 hot:
                # [row=1] ,,,,SKII15176
                # [row=2] ,,,文库描述,TB-QC06
                # [row=3] chr,start,end,mutation,SKII15177V1_1-TB-QC06-L04-AUDI27
                # 它的行 2 确实是文库描述，行 3 确实是文库全名。
                # 两个模式一致。
                for s in samples:
                    matrix[r_idx].append(desc_map.get(s, ""))
            elif r_idx == 2:
                for s in samples:
                    matrix[r_idx].append(s)
            else:
                row_vals = data_dict.get(key)
                vals = []
                if row_vals:
                    vals = row_vals[1:] if mode == "amp" else row_vals[4:]
                for s_idx, s in enumerate(samples):
                    if len(vals) > s_idx:
                        matrix[r_idx].append(vals[s_idx])
                    else:
                        matrix[r_idx].append("0")

    # 3. 整体清空，再整体写入
    clear_range(cfg["token"], cfg["sheet_name"], cfg["clean_range"])
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerows(matrix)
    put_csv_feishu(cfg["token"], cfg["sheet_name"], cfg["start_cell"], out.getvalue())
    
    # 4. 重新进行表头合并单元格以美化格式
    for start_let, end_let in merge_ranges:
        if start_let != end_let:
            print(f"[INFO] Merging header: {start_let}1:{end_let}1...")
            run_cmd(f"{LARK_CLI} sheets +cells-merge --spreadsheet-token '{cfg['token']}' --sheet-name '{cfg['sheet_name']}' --range '{start_let}1:{end_let}1'")

if __name__ == "__main__":
    main()
