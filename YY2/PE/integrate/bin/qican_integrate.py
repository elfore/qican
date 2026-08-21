#!/usr/bin/env python3
"""
One-command integration for YY2 PE batch summaries.

This tool keeps final summary CSV files as derived artifacts:
batch result directories -> full-history intermediate tables -> recomputed
concordance matrices -> project-specific summary CSV files.
"""

import argparse
import csv
import datetime as _dt
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path


BASE_DIR = Path("/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/PE")
INTEGRATE_DIR = BASE_DIR / "integrate"
CONFIG_DIR = INTEGRATE_DIR / "config"
MASTER_DIR = INTEGRATE_DIR / "master"
LOG_DIR = INTEGRATE_DIR / "logs"
WORK_DIR = INTEGRATE_DIR / "work"

DEFAULT_SUMMARY_DIR = Path(
    "/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/"
    "04.new_qican/S100/PE150/summary/YY2"
)
DEFAULT_BATCH_ROOT_PARENTS = [
    Path(
        "/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/"
        "04.new_qican/S100/PE150"
    ),
    Path(
        "/mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/"
        "04.new_qican/S100/PE150"
    ),
]
DEFAULT_SOURCE_INFO = BASE_DIR / "YY2_info.txt"
DEFAULT_KB_FILE = Path("/mnt/gpfs1/Users/caiyilun/test/qican_stat/chd/kownledge_base.txt")
DEFAULT_PB_POS = Path(
    "/mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/"
    "04.new_qican/S100/PE150/YY2-PL24-022-01/PB/result/combine_mut/"
    "PL2402201_LC002901-S01-D01-L12-C208C160_combine_rs.vcf"
)
DEFAULT_WGS_QC_PL = Path("/mnt/gpfs/Users/wangning/project/qican/WGS/pipeline/script/qc.pl")

PROJECTS = ("PA", "PB", "PD", "WGS")
DEFAULT_PROJECTS = ("PA", "PB", "PD")

BATCH_FIELDS = [
    "batch_id",
    "runID",
    "machineID",
    "version",
    "batch_root",
    "projects",
    "status",
]

SUMMARY_HEADERS = {
    ("PA", "01"): ["runID", "仪器号", "【其他信息】", "样本ID", "阳性符合率", "阴性符合率"],
    ("PB", "01"): [
        "runID",
        "仪器号",
        "【其他信息】",
        "样本ID",
        "阳性符合率",
        "阴性符合率",
        "CYP分型一致率",
        "HLA分型一致率",
        "HLA用药结果一致率",
    ],
    ("PD", "01"): [
        "runID",
        "仪器号",
        "【其他信息】",
        "样本ID",
        "阳性符合率",
        "阴性符合率",
        "CYP分型一致率",
        "HLA分型一致率",
        "HLA用药结果分型一致率",
    ],
    ("PA", "02"): [
        "runid",
        "mechineid",
        "libID",
        "raw_reads",
        "clean_reads",
        "effective_rate",
        "adapter",
        "raw_len",
        "filter_len",
        "Q20",
        "Q30",
        "GC",
        "mapped_reads",
        "mapped_rate",
        "coverage",
        "avg_depth",
        "on_target",
        "uniformity",
        "cov",
        "cov_20",
        "cov_100",
        "cov_200",
        "cov_500",
        "cov_1000",
        "var_num",
        "un_p_r",
        "cov",
        "cov_20",
        "cov_30",
        "cov_50",
        "cov_100",
        "cov_200",
        "cov_500",
        "cov_1000",
        "质控",
    ],
    ("PB", "02"): [
        "runinfo",
        "machineID",
        "libID",
        "x",
        "clean_reads",
        "effective_rate",
        "adapter",
        "raw_len",
        "filter_len",
        "Q20",
        "Q30",
        "GC",
        "mapped_reads",
        "mapped_rate",
        "coverage",
        "avg_depth",
        "on_target",
        "uniformity",
        "cov",
        "cov_10",
        "cov_20",
        "cov_30",
        "cov_40",
        "cov_50",
        "cov_100",
        "cov_200",
        "cov_500",
        "cov_1000",
        "warning_pos",
        "sample_name",
        "warning_info",
        "debug_info",
        "sample_judge",
    ],
    ("PD", "02"): [
        "runID",
        "mechineID",
        "sampleID",
        "raw_reads",
        "clean_reads",
        "effective_rate",
        "adapter",
        "raw_len",
        "filter_len",
        "Q20",
        "Q30",
        "GC",
        "mapped_reads",
        "mapped_rate",
        "coverage",
        "avg_depth",
        "on_target_whole",
        "on_target_hotspot",
        "on_target_cyp",
        "on_target_hla_a",
        "on_target_hla_b",
        "uniformity_snp",
        "uniformity_lg",
        "hotspot_cov",
        "hotspot_50_cov",
        "pb_cov",
        "pb_50_cov",
        "hla_a_cov",
        "hla_a_50_cov",
        "hla_b_cov",
        "hla_b_50_cov",
        "hla_c_cov",
        "hla_c_50_cov",
        "var_num",
        "un_p_r",
        "cov",
        "cov_20",
        "cov_30",
        "cov_50",
        "cov_100",
        "cov_200",
        "cov_500",
        "cov_1000",
        "是否合格",
    ],
}

SUMMARY_FILENAMES = {
    ("PA", "01"): "01_阴阳性符合率汇总.csv",
    ("PA", "02"): "02_质控汇总.csv",
    ("PA", "03"): "03_均一化覆盖度.csv",
    ("PA", "04"): "04_突变频率一致性.csv",
    ("PB", "01"): "01_阴阳性符合率汇总.csv",
    ("PB", "02"): "02_质控汇总.csv",
    ("PB", "03"): "03_均一化覆盖度.csv",
    ("PB", "04"): "04_突变频率一致性.csv",
    ("PB", "05"): "05_CYP分型结果汇总.csv",
    ("PB", "06"): "06_HLA分型结果汇总.csv",
    ("PB", "07"): "07_HLA用药结果汇总.csv",
    ("PD", "01"): "01_阴阳性符合率汇总.csv",
    ("PD", "02"): "02_质控汇总.csv",
    ("PD", "03"): "03_均一化覆盖度.csv",
    ("PD", "04"): "04_突变频率一致性.csv",
    ("PD", "05"): "05_CYP分型汇总.csv",
    ("PD", "06"): "06_HLA分型汇总.csv",
    ("PD", "07"): "07_HLA用药结果汇总.csv",
}

WGS_SPECIES_TO_FILE = {
    "Arabidopsis_thaliana": "01_拟南芥.csv",
    "Danio_rerio": "02_斑马鱼.csv",
    "Oryza_sativa": "03_水稻.csv",
    "Bos_taurus": "04_牛.csv",
}

DEFAULT_WGS_SAMPLE_SPECIES = {
    "CN002390": "Arabidopsis_thaliana",
    "CN002393": "Danio_rerio",
    "CN002389": "Oryza_sativa",
    "CN002392": "Bos_taurus",
}


class IntegrationError(RuntimeError):
    pass


def ensure_dirs():
    for path in (CONFIG_DIR, MASTER_DIR, LOG_DIR, WORK_DIR):
        path.mkdir(parents=True, exist_ok=True)
    for project in PROJECTS:
        (MASTER_DIR / project).mkdir(parents=True, exist_ok=True)


def now_tag():
    return _dt.datetime.now().strftime("%Y%m%d-%H%M%S")


def read_tsv(path, has_header=False):
    rows = []
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)
    if has_header and rows:
        return rows[0], rows[1:]
    return rows


def write_tsv(path, rows, header=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def write_csv(path, rows, header=None):
    rows = trim_empty_columns(rows, header=header)
    if header is not None:
        header, rows = rows[0], rows[1:]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def trim_empty_columns(rows, header=None):
    full_rows = []
    if header is not None:
        full_rows.append(list(header))
    full_rows.extend([list(row) for row in rows])
    if not full_rows:
        return []
    max_cols = max(len(row) for row in full_rows)
    keep_indexes = []
    for idx in range(max_cols):
        has_value = any(idx < len(row) and str(row[idx]).strip() != "" for row in full_rows)
        if has_value:
            keep_indexes.append(idx)
    if not keep_indexes:
        return [[] for _ in full_rows]
    trimmed = []
    for row in full_rows:
        trimmed.append([row[idx] if idx < len(row) else "" for idx in keep_indexes])
    return trimmed


def read_csv_rows(path):
    rows = []
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            rows.append(row)
    return rows


def csv_data_rows(path):
    rows = read_csv_rows(path)
    if not rows:
        return [], []
    return rows[0], [row for row in rows[1:] if row]


def history_rows_for_run(existing_csv, target_run_id):
    if not target_run_id or not existing_csv.exists():
        return []
    _, rows = csv_data_rows(existing_csv)
    return [row for row in rows if row_matches_run(row, target_run_id)]


def choose_target_rows(args, existing_csv, computed_rows):
    computed_target = [row for row in computed_rows if row_matches_run(row, args.run_id)]
    history_target = history_rows_for_run(existing_csv, args.run_id)
    source = getattr(args, "existing_target_source", "computed")
    if source == "history" and history_target:
        return history_target
    if source == "fallback" and not computed_target and history_target:
        return history_target
    return computed_target


def should_use_history_target(args, existing_csv, computed_rows):
    source = getattr(args, "existing_target_source", "computed")
    if source not in ("history", "fallback"):
        return False
    history_target = history_rows_for_run(existing_csv, args.run_id)
    if source == "history" and history_target:
        return True
    computed_target = [row for row in computed_rows if row_matches_run(row, args.run_id)]
    return source == "fallback" and not computed_target and bool(history_target)


def copy_if_exists(src, dst):
    if src.exists():
        try:
            if src.resolve() == dst.resolve():
                return
        except FileNotFoundError:
            pass
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(src), str(dst))


def copy_csv_clean_if_exists(src, dst):
    if src.exists():
        write_csv(dst, read_csv_rows(src))


def machine_from_run(run_id):
    parts = run_id.split("_")
    if len(parts) > 1 and parts[1].startswith("B"):
        return parts[1]
    return "-"


def batch_id_from_root(root):
    return Path(root).name


def batch_id_from_run_id(run_id):
    match = re.search(r"(SKII\d+)", run_id or "")
    if not match:
        raise IntegrationError("无法从 runID 推断 SKII 批次号: " + str(run_id))
    return match.group(1)


def infer_batch_root(run_id):
    batch_id = batch_id_from_run_id(run_id)
    candidates = [parent / batch_id for parent in DEFAULT_BATCH_ROOT_PARENTS]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def load_batches(path):
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in BATCH_FIELDS if field not in (reader.fieldnames or [])]
        if missing:
            raise IntegrationError("batches.tsv 缺少字段: " + ",".join(missing))
        return [dict(row) for row in reader if row.get("batch_root")]


def save_batches(path, batches):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BATCH_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in batches:
            writer.writerow({field: row.get(field, "") for field in BATCH_FIELDS})


def init_batches_from_source(source_info):
    rows = []
    with open(source_info, "r", encoding="utf-8-sig") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            run_id, version, root = parts[:3]
            root_path = Path(root)
            projects = [p for p in PROJECTS if (root_path / p).exists()]
            if not projects:
                projects = ["PA", "PB", "PD"]
            rows.append(
                {
                    "batch_id": batch_id_from_root(root),
                    "runID": run_id,
                    "machineID": machine_from_run(run_id),
                    "version": version,
                    "batch_root": root,
                    "projects": ",".join(projects),
                    "status": "done",
                }
            )
    return rows


def merge_batch(batches, new_row):
    key = new_row["batch_root"]
    replaced = False
    out = []
    for row in batches:
        if row["batch_root"] == key or row["runID"] == new_row["runID"]:
            out.append(new_row)
            replaced = True
        else:
            out.append(row)
    if not replaced:
        out.append(new_row)
    return out


def active_batches(batches, projects):
    wanted = set(projects)
    out = []
    for row in batches:
        if row.get("status", "done") != "done":
            continue
        row_projects = set(x.strip() for x in row.get("projects", "").split(",") if x.strip())
        if row_projects & wanted:
            out.append(row)
    return out


def write_info_file(path, batches, projects=None):
    rows = []
    wanted = set(projects or PROJECTS)
    for row in batches:
        row_projects = set(x.strip() for x in row.get("projects", "").split(",") if x.strip())
        if row.get("status", "done") != "done":
            continue
        if projects and not (row_projects & wanted):
            continue
        rows.append([row["runID"], row.get("version") or row.get("machineID") or "-", row["batch_root"]])
    write_tsv(path, rows)
    return len(rows)


def select_dp_batches(args, batches):
    if args.dp_scope == "skip":
        return []
    scoped = [row for row in batches if "_" in row.get("runID", "")]
    if args.dp_scope == "target" and getattr(args, "run_id", ""):
        scoped = [row for row in scoped if row.get("runID") == args.run_id]
    return scoped


def target_batches(args, batches):
    run_id = getattr(args, "run_id", "")
    if not run_id:
        return batches
    return [row for row in batches if row.get("runID") == run_id]


def stats_batches(args, batches):
    if getattr(args, "history_mode", "summary") != "summary" or not getattr(args, "run_id", ""):
        return batches
    selected = []
    seen = set()
    for row in batches:
        if row.get("runID") in ("PL2402201", args.run_id) and row.get("runID") not in seen:
            selected.append(row)
            seen.add(row.get("runID"))
    return selected


def history_summary_path(args, project, key):
    return Path(args.history_summary_dir) / project / SUMMARY_FILENAMES[(project, key)]


def run_cmd(cmd, cwd=None, dry_run=False):
    if dry_run:
        print("[dry-run:exec-temp] " + " ".join(map(str, cmd)), flush=True)
        return
    else:
        print("[run] " + " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], cwd=str(cwd) if cwd else None, check=True)


def first_existing(patterns):
    for pattern in patterns:
        matches = sorted(Path().glob(pattern) if not pattern.startswith("/") else Path("/").glob(pattern[1:]))
        if matches:
            return matches[0]
    return None


def glob_paths(pattern):
    if pattern.startswith("/"):
        return sorted(Path("/").glob(pattern[1:]))
    return sorted(Path().glob(pattern))


def find_pd_pos(batch_root):
    matches = glob_paths(str(Path(batch_root) / "PD" / "result" / "mutation" / "*final_mut.txt"))
    return matches[0] if matches else None


def project_work_dir(run_work, project):
    path = run_work / project
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_existing_project_scripts(args, batches, run_work):
    info_file = run_work / "YY2_info.full.tsv"
    dp_info_file = run_work / "YY2_info.dp.tsv"
    calc_batches = stats_batches(args, batches)
    write_info_file(info_file, calc_batches, projects=["PA", "PB", "PD"])
    dp_count = write_info_file(dp_info_file, select_dp_batches(args, batches), projects=["PA", "PB", "PD"])
    py = sys.executable

    if "PA" in args.projects:
        outdir = project_work_dir(run_work, "PA")
        bin_dir = BASE_DIR / "PA" / "bin"
        run_cmd([py, bin_dir / "stat.py", "-infile", info_file, "-pos", args.pa_pos, "-outfile", "qc_stat.txt", "--outdir", outdir], dry_run=args.dry_run)
        run_cmd([py, bin_dir / "sample_stat.py", "-infile", info_file, "-pos", args.pa_pos, "-outfile", "vars_stat.txt", "--outdir", outdir], dry_run=args.dry_run)
        run_cmd([py, bin_dir / "sample_stat_allv.py", "-infile", info_file, "-pos", args.pa_pos, "-outfile", "vars_stat_all.txt", "--outdir", outdir], dry_run=args.dry_run)
        if dp_count:
            run_cmd([py, bin_dir / "stat_dp_1.py", "-infile", dp_info_file, "-pos", args.pa_pos, "-outfile", "dp_stat.txt", "--outdir", outdir], dry_run=args.dry_run)

    if "PB" in args.projects:
        outdir = project_work_dir(run_work, "PB")
        bin_dir = BASE_DIR / "PB" / "bin"
        hla_file = bin_dir / "HLA.txt"
        run_cmd([py, bin_dir / "stat.py", "-infile", info_file, "-outfile", "result.txt", "--outdir", outdir], dry_run=args.dry_run)
        run_cmd([py, bin_dir / "sample_stat.py", "-infile", info_file, "-pos", args.pb_pos, "-outfile", "vars_stat.txt", "--outdir", outdir], dry_run=args.dry_run)
        run_cmd([py, bin_dir / "sample_stat_allv.py", "-infile", info_file, "-pos", args.pb_pos, "-outfile", "vars_stat_all.txt", "-infile2", hla_file, "--outdir", outdir], dry_run=args.dry_run)
        if dp_count:
            run_cmd([py, bin_dir / "stat_dp_1.py", "-infile", dp_info_file, "-pos", args.pb_pos, "-outfile", "dp_stat.txt", "--outdir", outdir], dry_run=args.dry_run)

    if "PD" in args.projects:
        outdir = project_work_dir(run_work, "PD")
        bin_dir = BASE_DIR / "PD" / "bin"
        hla_file = bin_dir / "HLA.txt"
        pd_pos = args.pd_pos or find_pd_pos(args.batch_root)
        if not pd_pos:
            print("[warn] 未找到 PD -pos 文件，跳过 PD vars/dp 相关脚本")
        else:
            run_cmd([py, bin_dir / "stat.py", "-infile", info_file, "-pos", pd_pos, "-outfile", "qc_stat.txt", "--outdir", outdir], dry_run=args.dry_run)
            run_cmd([py, bin_dir / "sample_stat.py", "-infile", info_file, "-pos", pd_pos, "-outfile", "vars_stat.txt", "--outdir", outdir], dry_run=args.dry_run)
            run_cmd([py, bin_dir / "sample_stat_allv.py", "-infile", info_file, "-pos", pd_pos, "-outfile", "vars_stat_all.txt", "-infile2", hla_file, "--outdir", outdir], dry_run=args.dry_run)
            if dp_count:
                run_cmd([py, bin_dir / "stat_dp_1.py", "-infile", dp_info_file, "-pos", pd_pos, "-outfile", "dp_stat.txt", "--outdir", outdir], dry_run=args.dry_run)


def normalize_sample_from_cyp(path, run_id):
    name = Path(path).name.split(".")[0]
    sample_name2 = name
    sample_name = name.split("-")[0]
    if "_" in sample_name2:
        if run_id == "PL2402201":
            sample_name2 = sample_name2.split("_", 1)[1]
            sample_name = sample_name2.split("-")[0]
        else:
            rhs = sample_name2.split("_", 1)[1]
            sample_name2 = "-".join(rhs.split("-")[1:])
            sample_name = sample_name2.split("-")[0]
    return sample_name, sample_name2


def normalize_sample_from_hla(path, run_id):
    stem = Path(path).name.split("_hla")[0]
    sample_name2 = stem
    sample_name = stem.split("-")[0]
    if "_" in sample_name2:
        if run_id == "PL2402201":
            sample_name2 = sample_name2.split("_", 1)[1]
            sample_name = sample_name2.split("-")[0]
        else:
            rhs = sample_name2.split("_", 1)[1]
            sample_name2 = "-".join(rhs.split("-")[1:])
            sample_name = sample_name2.split("-")[0]
    return sample_name, sample_name2


def parse_hla_drug_map(hla_file):
    maps = {"A": {}, "B": {}, "C": {}}
    if not hla_file.exists():
        return maps
    with open(hla_file, "r", encoding="utf-8-sig") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            gene = parts[3].replace("HLA-", "")
            if gene in maps:
                maps[gene].setdefault(parts[0], parts[1] + "|" + parts[5])
    return maps


def collect_cyp_typing(project, batches, target_run_id=None):
    data = defaultdict(lambda: defaultdict(dict))
    order = defaultdict(list)
    for batch in batches:
        if project not in batch.get("projects", ""):
            continue
        root = Path(batch["batch_root"]) / project / "result"
        if project == "PB":
            patterns = [root / "cyp_analysis" / "*cyp_allele.txt", root / "cyp_analysis_new" / "*cyp_allele.txt"]
        else:
            patterns = [root / "genotyping" / "*allele.txt"]
        files = []
        for pattern in patterns:
            files.extend(sorted(pattern.parent.glob(pattern.name)))
        for file_path in files:
            sample, sample2 = normalize_sample_from_cyp(file_path, batch["runID"])
            test_id = batch["runID"] + "_" + sample2
            if test_id not in order[sample]:
                order[sample].append(test_id)
            with open(file_path, "r", encoding="utf-8-sig") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    data[sample][parts[0]][test_id] = parts[1]
    return matrix_majority(data, order, target_run_id)


def collect_hla_typing(project, batches, target_run_id=None):
    data = defaultdict(lambda: defaultdict(dict))
    order = defaultdict(list)
    for batch in batches:
        if project not in batch.get("projects", ""):
            continue
        result_dir = Path(batch["batch_root"]) / project / "result" / "hla_analysis_new2"
        for file_path in sorted(result_dir.glob("*hla_result.tsv")):
            sample, sample2 = normalize_sample_from_hla(file_path, batch["runID"])
            test_id = batch["runID"] + "_" + sample2
            if test_id not in order[sample]:
                order[sample].append(test_id)
            with open(file_path, "r", encoding="utf-8-sig") as handle:
                for line in handle:
                    if line.startswith("Gene"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    key = parts[0].split("_")[0][-1] + parts[0][-1]
                    data[sample][key][test_id] = parts[1]
    align_hla_pairs(data, order)
    return matrix_baseline(data, order, target_run_id)


def collect_hla_drug(project, batches, hla_file, target_run_id=None):
    drug_map = parse_hla_drug_map(hla_file)
    data = defaultdict(lambda: defaultdict(dict))
    order = defaultdict(list)
    for batch in batches:
        if project not in batch.get("projects", ""):
            continue
        result_dir = Path(batch["batch_root"]) / project / "result" / "hla_analysis_new2"
        for file_path in sorted(result_dir.glob("*hla_result.tsv")):
            sample, sample2 = normalize_sample_from_hla(file_path, batch["runID"])
            if sample2 in ("CN001799-S08-D01-L09-C220C159", "CN001799-S08-D01-L10-C220C163"):
                continue
            test_id = batch["runID"] + "_" + sample2
            if test_id not in order[sample]:
                order[sample].append(test_id)
            tmp = defaultdict(list)
            with open(file_path, "r", encoding="utf-8-sig") as handle:
                for line in handle:
                    if line.startswith("Gene"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    gene = parts[0].split("_")[0][-1]
                    allele = ":".join(parts[1].split(":")[0:2])[1:]
                    if gene in ("A", "B", "C"):
                        tmp[gene].append(allele)
            for gene, alleles in tmp.items():
                values = [drug_map[gene][allele] for allele in alleles if allele in drug_map.get(gene, {})]
                value = "阴性" if not values else "----".join(sorted(set(values)))
                data[sample]["HLA-" + gene][test_id] = value
    return matrix_majority(data, order, target_run_id)


def target_samples(order, target_run_id):
    if not target_run_id:
        return set(order)
    return {sample for sample, ids in order.items() if any(x.startswith(target_run_id + "_") for x in ids)}


def matrix_majority(data, order, target_run_id=None):
    rows = []
    keep = target_samples(order, target_run_id)
    for sample in sorted(data):
        if sample not in keep:
            continue
        ids = order[sample]
        if not ids:
            continue
        rows.append([sample, "一致率"] + ids)
        for key in sorted(data[sample]):
            values = [data[sample][key].get(test_id, "NA") for test_id in ids]
            counts = Counter(values)
            rate = "NA" if not values else "%.2f%%" % (max(counts.values()) / len(values) * 100)
            rows.append([key, rate] + values)
        rows.append([])
    return rows


def align_hla_pairs(data, order):
    for sample, sample_data in data.items():
        ids = order[sample]
        if not ids:
            continue
        baseline = ids[0]
        for gene in sorted({key[:-1] for key in sample_data if len(key) > 1}):
            v1, v2 = gene + "1", gene + "2"
            if v1 not in sample_data or v2 not in sample_data:
                continue
            b1 = sample_data[v1].get(baseline, "NA")
            b2 = sample_data[v2].get(baseline, "NA")
            for test_id in ids[1:]:
                t1 = sample_data[v1].get(test_id, "NA")
                t2 = sample_data[v2].get(test_id, "NA")
                score_no_swap = int(t1 == b1 and b1 != "NA") + int(t2 == b2 and b2 != "NA")
                score_swap = int(t2 == b1 and b1 != "NA") + int(t1 == b2 and b2 != "NA")
                if score_swap > score_no_swap:
                    sample_data[v1][test_id] = t2
                    sample_data[v2][test_id] = t1


def matrix_baseline(data, order, target_run_id=None):
    rows = []
    keep = target_samples(order, target_run_id)
    for sample in sorted(data):
        if sample not in keep:
            continue
        ids = order[sample]
        if not ids:
            continue
        baseline = ids[0]
        tests = ids[1:]
        rows.append([sample, "一致率"] + ids)
        for key in sorted(data[sample]):
            values = [data[sample][key].get(test_id, "NA") for test_id in ids]
            b_val = data[sample][key].get(baseline, "NA")
            if tests:
                matches = sum(1 for test_id in tests if data[sample][key].get(test_id, "NA") == b_val and b_val != "NA")
                rate = "%d%%" % int(matches / len(tests) * 100)
            else:
                rate = "NA"
            rows.append([key, rate] + values)
        rows.append([])
    return rows


def render_01(project, source_tsv, out_csv):
    rendered = render_01_rows(project, source_tsv)
    write_csv(out_csv, rendered, header=SUMMARY_HEADERS[(project, "01")])
    return len(rendered)


def render_01_rows(project, source_tsv):
    rows = read_tsv(source_tsv)
    if not rows:
        return []
    data_rows = rows[1:] if rows[0] and rows[0][0] == "runID" else rows
    rendered = []
    for row in data_rows:
        if not row:
            continue
        base = [
            row[0] if len(row) > 0 else "",
            row[1] if len(row) > 1 else "",
            "",
            row[2] if len(row) > 2 else "",
            row[3] if len(row) > 3 else "",
            row[4] if len(row) > 4 else "",
        ]
        if project in ("PB", "PD"):
            if len(row) > 7:
                base.extend([
                    row[5] if len(row) > 5 else "",
                    row[6] if len(row) > 6 else "",
                    row[7] if len(row) > 7 else "",
                ])
            else:
                # Legacy PB/PD sample_stat_allv.py writes only HLA typing and HLA-drug
                # metrics after positive/negative rates. CYP is filled from 05 matrix.
                base.extend([
                    "",
                    row[5] if len(row) > 5 else "",
                    row[6] if len(row) > 6 else "",
                ])
        rendered.append(base)
    return rendered


def render_tsv_with_header(project, key, source_tsv, out_csv):
    header, rows = read_tsv(source_tsv, has_header=True)
    write_csv(out_csv, rows, header=header)
    return len(rows)


def tsv_rows_with_header(source_tsv):
    header, rows = read_tsv(source_tsv, has_header=True)
    return header, [row for row in rows if row]


def tsv_rows_without_header(source_tsv):
    return [row for row in read_tsv(source_tsv) if row]


def row_matches_run(row, run_id):
    return bool(run_id) and bool(row) and row[0] == run_id


def blank_summary_key(row, key_indexes):
    return not any((row[i].strip() if i < len(row) else "") for i in key_indexes)


def merge_summary_rows(existing_csv, out_csv, header, new_rows, key_indexes, target_run_id=None):
    merged = OrderedDict()
    if existing_csv.exists():
        old_header, old_rows = csv_data_rows(existing_csv)
        if old_header:
            header = old_header
        for row in old_rows:
            if blank_summary_key(row, key_indexes):
                continue
            key = tuple(row[i] if i < len(row) else "" for i in key_indexes)
            merged[key] = row
    kept_new = 0
    for row in new_rows:
        if target_run_id and not row_matches_run(row, target_run_id):
            continue
        if blank_summary_key(row, key_indexes):
            continue
        key = tuple(row[i] if i < len(row) else "" for i in key_indexes)
        merged[key] = row
        kept_new += 1
    write_csv(out_csv, list(merged.values()), header=header)
    return kept_new


def merge_tsv_rows_with_existing(existing_csv, out_csv, header, new_rows, key_indexes=(0, 2)):
    merged = OrderedDict()
    if existing_csv.exists():
        old_rows = read_csv_rows(existing_csv)
        if old_rows:
            if old_rows[0]:
                header = old_rows[0]
            for row in old_rows[1:]:
                if not row or blank_summary_key(row, key_indexes):
                    continue
                key = tuple(row[i] if i < len(row) else "" for i in key_indexes)
                merged[key] = row
    for row in new_rows:
        if not row or blank_summary_key(row, key_indexes):
            continue
        key = tuple(row[i] if i < len(row) else "" for i in key_indexes)
        merged[key] = row
    write_csv(out_csv, list(merged.values()), header=header)
    return len(new_rows)


def render_tsv_merge_existing(source_tsv, existing_csv, out_csv, key_indexes=(0, 2)):
    header, new_rows = read_tsv(source_tsv, has_header=True)
    return merge_tsv_rows_with_existing(existing_csv, out_csv, header, [row for row in new_rows if row], key_indexes)


def count_tsv_data_rows(source_tsv):
    count = 0
    with open(source_tsv, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for idx, row in enumerate(reader):
            if idx == 0:
                continue
            if row:
                count += 1
    return count


def render_headerless(project, key, source_tsv, out_csv):
    rows = read_tsv(source_tsv)
    write_csv(out_csv, rows, header=SUMMARY_HEADERS[(project, key)])
    return len(rows)


def trivial_na_matrix_item(row, max_values=2):
    if not row:
        return False
    if row[0].strip():
        values = [value.strip() for value in row[2:] if value.strip()]
    else:
        values = [value.strip() for value in row if value.strip()]
    return bool(values) and len(values) <= max_values and all(value == "NA" for value in values)


def trivial_rate_matrix_item(row):
    if not row or row[0].strip():
        return False
    rate = row[1].strip() if len(row) > 1 else ""
    tail = [value.strip() for value in row[2:] if value.strip()]
    return rate in ("100%", "100.00%") and all(value == "NA" for value in tail)


def filter_matrix_rows(rows, drop_trivial_na_items=False, drop_trivial_rate_items=False):
    if not drop_trivial_na_items and not drop_trivial_rate_items:
        return rows
    filtered = []
    for row in rows:
        if drop_trivial_na_items and trivial_na_matrix_item(row):
            continue
        if drop_trivial_rate_items and trivial_rate_matrix_item(row):
            continue
        filtered.append(row)
    return filtered


def write_matrix(project, key, rows, out_csv, master_tsv):
    rows = filter_matrix_rows(
        rows,
        drop_trivial_na_items=(key == "06"),
        drop_trivial_rate_items=(key in ("05", "07")),
    )
    write_csv(out_csv, rows)
    write_tsv(master_tsv, rows)
    return sum(1 for row in rows if row)


def format_fraction(value):
    return str(value)


def matrix_baseline_scores_for_run(matrix_csv, target_run_id):
    rows = read_csv_rows(matrix_csv)
    scores = {}
    i = 0
    while i < len(rows):
        row = rows[i]
        if len(row) >= 2 and row[1] == "一致率":
            header = row
            block = []
            i += 1
            while i < len(rows) and rows[i]:
                if rows[i][0].strip():
                    block.append(rows[i])
                i += 1
            for idx, test_id in enumerate(header):
                if idx < 2 or not test_id.startswith(target_run_id + "_"):
                    continue
                matches = 0
                total = 0
                for item in block:
                    baseline = item[2] if len(item) > 2 else "NA"
                    value = item[idx] if idx < len(item) else "NA"
                    if value == "":
                        continue
                    total += 1
                    if value == baseline:
                        matches += 1
                sample_id = test_id[len(target_run_id) + 1:]
                scores[sample_id] = format_fraction(matches / total) if total else ""
        i += 1
    return scores


def update_01_cyp_from_matrix(project, run_id, matrix_csv, summary_csv):
    if project not in ("PB", "PD") or not matrix_csv.exists() or not summary_csv.exists():
        return 0
    scores = matrix_baseline_scores_for_run(matrix_csv, run_id)
    if not scores:
        return 0
    rows = read_csv_rows(summary_csv)
    updated = 0
    for row in rows[1:]:
        if not row or row[0] != run_id:
            continue
        while len(row) < 9:
            row.append("")
        sample_id = row[3] if len(row) > 3 else ""
        if sample_id in scores:
            row[6] = scores[sample_id]
            updated += 1
    write_csv(summary_csv, rows)
    return updated


def split_matrix_blocks(rows):
    blocks = OrderedDict()
    current_sample = None
    for row in rows:
        if not row:
            current_sample = None
            continue
        if len(row) >= 2 and row[1] == "一致率":
            current_sample = row[0]
            blocks[current_sample] = {"header": row[:], "items": OrderedDict()}
            continue
        if current_sample and row[0]:
            blocks[current_sample]["items"][row[0]] = row[:]
    return blocks


def matrix_rows_from_blocks(blocks, drop_trivial_na_items=False, drop_trivial_rate_items=False):
    rows = []
    for block in blocks.values():
        rows.append(block["header"])
        for item in block["items"].values():
            if drop_trivial_na_items and trivial_na_matrix_item(item):
                continue
            if drop_trivial_rate_items and trivial_rate_matrix_item(item):
                continue
            rows.append(item)
        rows.append([])
    return rows


def recompute_matrix_rate(values, mode):
    clean = [value for value in values if value != ""]
    if not clean:
        return "NA"
    if mode == "baseline":
        baseline = clean[0]
        tests = clean[1:]
        if not tests:
            return "NA"
        matches = sum(1 for value in tests if value == baseline and baseline != "NA")
        return "%d%%" % int(matches / len(tests) * 100)
    counts = Counter(clean)
    return "%.2f%%" % (max(counts.values()) / len(clean) * 100)


def align_hla_block_to_history(block, ids, target_ids):
    if not ids or not target_ids:
        return
    id_index = {test_id: idx for idx, test_id in enumerate(ids)}
    target_indexes = [id_index[test_id] for test_id in target_ids if test_id in id_index]
    if not target_indexes:
        return
    item_keys = set(block["items"].keys())
    for gene in sorted({key[:-1] for key in item_keys if len(key) > 1}):
        v1, v2 = gene + "1", gene + "2"
        if v1 not in block["items"] or v2 not in block["items"]:
            continue
        row1 = block["items"][v1]
        row2 = block["items"][v2]
        b1 = row1[2] if len(row1) > 2 else "NA"
        b2 = row2[2] if len(row2) > 2 else "NA"
        for idx in target_indexes:
            pos = idx + 2
            t1 = row1[pos] if pos < len(row1) else "NA"
            t2 = row2[pos] if pos < len(row2) else "NA"
            score_no_swap = int(t1 == b1 and b1 != "NA") + int(t2 == b2 and b2 != "NA")
            score_swap = int(t2 == b1 and b1 != "NA") + int(t1 == b2 and b2 != "NA")
            if score_swap > score_no_swap:
                row1[pos], row2[pos] = t2, t1


def refresh_matrix_rates(block, mode):
    for item_key, row in block["items"].items():
        values = row[2:]
        row[1] = recompute_matrix_rate(values, mode)


def matrix_contains_run(existing_csv, target_run_id):
    if not target_run_id or not existing_csv.exists():
        return False
    for row in read_csv_rows(existing_csv):
        if len(row) >= 3 and row[1] == "一致率":
            if any(test_id.startswith(target_run_id + "_") for test_id in row[2:]):
                return True
    return False


def merge_matrix_with_history(existing_csv, new_rows, out_csv, mode, drop_trivial_na_items=False, drop_trivial_rate_items=False):
    history_rows = read_csv_rows(existing_csv) if existing_csv.exists() else []
    blocks = split_matrix_blocks(history_rows)
    new_blocks = split_matrix_blocks(new_rows)
    updated = 0
    for sample, new_block in new_blocks.items():
        if sample not in blocks:
            blocks[sample] = new_block
            updated += sum(1 for row in new_block["items"].values() if row)
            continue
        block = blocks[sample]
        old_ids = block["header"][2:]
        new_ids = new_block["header"][2:]
        ids = old_ids[:]
        for test_id in new_ids:
            if test_id not in ids:
                ids.append(test_id)
        block["header"] = [sample, "一致率"] + ids
        item_keys = list(block["items"].keys())
        for item_key in new_block["items"]:
            if item_key not in item_keys:
                item_keys.append(item_key)
        for item_key in item_keys:
            old_values = {}
            if item_key in block["items"]:
                old_row = block["items"][item_key]
                for idx, test_id in enumerate(old_ids):
                    old_values[test_id] = old_row[idx + 2] if idx + 2 < len(old_row) else "NA"
            if item_key in new_block["items"]:
                new_item = new_block["items"][item_key]
                for idx, test_id in enumerate(new_ids):
                    old_values[test_id] = new_item[idx + 2] if idx + 2 < len(new_item) else "NA"
            else:
                for test_id in new_ids:
                    old_values[test_id] = "NA"
            values = [old_values.get(test_id, "NA") for test_id in ids]
            rate = recompute_matrix_rate(values, mode)
            block["items"][item_key] = [item_key, rate] + values
            updated += 1
        if mode == "baseline":
            align_hla_block_to_history(block, ids, new_ids)
            refresh_matrix_rates(block, mode)
    rows = matrix_rows_from_blocks(
        blocks,
        drop_trivial_na_items=drop_trivial_na_items,
        drop_trivial_rate_items=drop_trivial_rate_items,
    )
    write_csv(out_csv, rows)
    return sum(1 for row in rows if row)


def should_write_preview(args):
    return bool(args.write or args.preview_files)


def render_project_outputs(args, batches, run_work):
    summary_dir = Path(args.summary_dir)
    counts = OrderedDict()
    if args.write and not args.dry_run:
        backup_summary(summary_dir)

    for project in ("PA", "PB", "PD"):
        if project not in args.projects:
            continue
        project_out = summary_dir / project
        project_out.mkdir(parents=True, exist_ok=True)
        master_out = MASTER_DIR / project
        work = run_work / project

        source = work / "vars_stat_all.txt"
        if source.exists():
            dst = project_out / SUMMARY_FILENAMES[(project, "01")]
            print("[render] %s 01 <- %s" % (project, source), flush=True)
            if should_write_preview(args):
                target = dst if args.write else run_work / "preview" / project / dst.name
                existing = history_summary_path(args, project, "01")
                if args.history_mode == "summary":
                    computed_rows = render_01_rows(project, source)
                    if should_use_history_target(args, existing, computed_rows):
                        copy_csv_clean_if_exists(existing, target)
                        counts[(project, "01")] = len(history_rows_for_run(existing, args.run_id))
                    else:
                        rows = choose_target_rows(args, existing, computed_rows)
                        counts[(project, "01")] = merge_summary_rows(
                            existing,
                            target,
                            SUMMARY_HEADERS[(project, "01")],
                            rows,
                            key_indexes=(0, 3),
                        )
                else:
                    counts[(project, "01")] = render_01(project, source, target)
            else:
                counts[(project, "01")] = count_tsv_data_rows(source)
            if args.write:
                copy_if_exists(source, master_out / "positive_negative_source.tsv")

        source = work / ("result.txt" if project == "PB" else "qc_stat.txt")
        if source.exists():
            dst = project_out / SUMMARY_FILENAMES[(project, "02")]
            print("[render] %s 02 <- %s" % (project, source), flush=True)
            if should_write_preview(args):
                target = dst if args.write else run_work / "preview" / project / dst.name
                existing = history_summary_path(args, project, "02")
                if args.history_mode == "summary":
                    computed_rows = tsv_rows_without_header(source)
                    if should_use_history_target(args, existing, computed_rows):
                        copy_csv_clean_if_exists(existing, target)
                        counts[(project, "02")] = len(history_rows_for_run(existing, args.run_id))
                    else:
                        rows = choose_target_rows(args, existing, computed_rows)
                        counts[(project, "02")] = merge_summary_rows(
                            existing,
                            target,
                            SUMMARY_HEADERS[(project, "02")],
                            rows,
                            key_indexes=(0, 2),
                        )
                else:
                    counts[(project, "02")] = render_headerless(project, "02", source, target)
            else:
                counts[(project, "02")] = len(read_tsv(source))
            if args.write:
                copy_if_exists(source, master_out / "qc.tsv")

        source = work / "dp_stat.txt"
        if source.exists():
            dst = project_out / SUMMARY_FILENAMES[(project, "03")]
            print("[render] %s 03 <- %s" % (project, source), flush=True)
            target = dst if args.write else run_work / "preview" / project / dst.name
            if not should_write_preview(args):
                counts[(project, "03")] = count_tsv_data_rows(source)
                print("[dry-run] %s 03 counted only, skip writing wide coverage preview" % project, flush=True)
            elif args.dp_scope == "target" and not args.write:
                counts[(project, "03")] = count_tsv_data_rows(source)
                print("[dry-run] %s 03 counted only, skip writing wide coverage preview" % project, flush=True)
            elif args.dp_scope == "target" and args.write:
                header, rows = tsv_rows_with_header(source)
                existing = history_summary_path(args, project, "03")
                if args.history_mode == "summary" and should_use_history_target(args, existing, rows):
                    copy_csv_clean_if_exists(existing, target)
                    counts[(project, "03")] = len(history_rows_for_run(existing, args.run_id))
                else:
                    rows = choose_target_rows(args, existing, rows) if args.history_mode == "summary" else rows
                    counts[(project, "03")] = merge_tsv_rows_with_existing(existing, target, header, rows)
            else:
                counts[(project, "03")] = render_tsv_with_header(project, "03", source, target)
            if args.write:
                copy_if_exists(source, master_out / "depth.tsv")

        source = work / "vars_stat.txt"
        if source.exists():
            dst = project_out / SUMMARY_FILENAMES[(project, "04")]
            print("[render] %s 04 <- %s" % (project, source), flush=True)
            if should_write_preview(args):
                target = dst if args.write else run_work / "preview" / project / dst.name
                existing = history_summary_path(args, project, "04")
                if args.history_mode == "summary":
                    header, rows = tsv_rows_with_header(source)
                    if should_use_history_target(args, existing, rows):
                        copy_csv_clean_if_exists(existing, target)
                        counts[(project, "04")] = len(history_rows_for_run(existing, args.run_id))
                    else:
                        rows = choose_target_rows(args, existing, rows)
                        counts[(project, "04")] = merge_summary_rows(
                            existing,
                            target,
                            header,
                            rows,
                            key_indexes=(0, 2),
                        )
                else:
                    counts[(project, "04")] = render_tsv_with_header(project, "04", source, target)
            else:
                counts[(project, "04")] = count_tsv_data_rows(source)
            if args.write:
                copy_if_exists(source, master_out / "var_freq.tsv")

    for project in ("PB", "PD"):
        if project not in args.projects:
            continue
        project_out = summary_dir / project
        project_out.mkdir(parents=True, exist_ok=True)
        master_out = MASTER_DIR / project
        hla_file = BASE_DIR / project / "bin" / "HLA.txt"
        if args.history_mode == "summary":
            print("[concordance] %s CYP/HLA/HLA-drug from summary history + target result files" % project, flush=True)
            concordance_batches = target_batches(args, batches)
        else:
            print("[concordance] %s CYP/HLA/HLA-drug from result files" % project, flush=True)
            concordance_batches = batches
        matrices = {
            "05": collect_cyp_typing(project, concordance_batches, target_run_id=args.run_id),
            "06": collect_hla_typing(project, concordance_batches, target_run_id=args.run_id),
            "07": collect_hla_drug(project, concordance_batches, hla_file, target_run_id=args.run_id),
        }
        master_names = {"05": "cyp_typing.tsv", "06": "hla_typing.tsv", "07": "hla_drug.tsv"}
        matrix_modes = {"05": "majority", "06": "baseline", "07": "majority"}
        for key, rows in matrices.items():
            dst = project_out / SUMMARY_FILENAMES[(project, key)]
            preview_dst = run_work / "preview" / project / dst.name
            if should_write_preview(args):
                target = dst if args.write else preview_dst
                existing = history_summary_path(args, project, key)
                if args.history_mode == "summary":
                    use_history_target = (
                        getattr(args, "existing_target_source", "computed") == "history"
                        and matrix_contains_run(existing, args.run_id)
                    )
                    fallback_history_target = (
                        getattr(args, "existing_target_source", "computed") == "fallback"
                        and not rows
                        and matrix_contains_run(existing, args.run_id)
                    )
                    if use_history_target or fallback_history_target:
                        copy_csv_clean_if_exists(existing, target)
                        counts[(project, key)] = sum(1 for row in read_csv_rows(existing) if row)
                    else:
                        counts[(project, key)] = merge_matrix_with_history(
                            existing,
                            rows,
                            target,
                            matrix_modes[key],
                            drop_trivial_na_items=(key == "06"),
                            drop_trivial_rate_items=(key in ("05", "07")),
                        )
                    if args.write:
                        copy_if_exists(target, master_out / master_names[key])
                else:
                    counts[(project, key)] = write_matrix(
                        project,
                        key,
                        rows,
                        target,
                        master_out / master_names[key] if args.write else run_work / "preview" / project / master_names[key],
                    )
            else:
                counts[(project, key)] = sum(1 for row in rows if row)
        if should_write_preview(args):
            matrix05 = project_out / SUMMARY_FILENAMES[(project, "05")] if args.write else run_work / "preview" / project / SUMMARY_FILENAMES[(project, "05")]
            summary01 = project_out / SUMMARY_FILENAMES[(project, "01")] if args.write else run_work / "preview" / project / SUMMARY_FILENAMES[(project, "01")]
            fixed = update_01_cyp_from_matrix(project, args.run_id, matrix05, summary01)
            if fixed:
                print("[render] %s 01 CYP <- %s (%s rows)" % (project, matrix05, fixed), flush=True)

    return counts


def backup_summary(summary_dir):
    if not summary_dir.exists():
        return
    backup_dir = summary_dir.parent / ("YY2.backup." + now_tag())
    print("[backup] " + str(summary_dir) + " -> " + str(backup_dir))
    shutil.copytree(str(summary_dir), str(backup_dir))


def check_batch(args, batches):
    target = Path(args.batch_root)
    issues = []
    if not target.exists():
        issues.append("批次目录不存在: " + str(target))
    for project in args.projects:
        if project == "WGS":
            continue
        pdir = target / project
        if not pdir.exists():
            issues.append(project + " 目录不存在: " + str(pdir))
    if "PB" in args.projects and (target / "PB").exists() and not (target / "PB" / "result").exists():
        issues.append("PB/result 不存在，可能仍停留在 Nextflow work 目录，PB 汇总可能为空")
    if "PD" in args.projects and not args.pd_pos and not find_pd_pos(target):
        issues.append("PD 未找到 result/mutation/*final_mut.txt，PD 阴阳性/突变频率相关统计可能无法生成")
    return issues


def parse_projects(value):
    projects = [x.strip() for x in value.split(",") if x.strip()]
    bad = [x for x in projects if x not in PROJECTS]
    if bad:
        raise argparse.ArgumentTypeError("未知项目: " + ",".join(bad))
    return projects


def command_init_batches(args):
    ensure_dirs()
    rows = init_batches_from_source(Path(args.source_info))
    save_batches(Path(args.batches), rows)
    print("已初始化批次清单: %s, %d rows" % (args.batches, len(rows)))


def command_add_batch(args):
    ensure_dirs()
    if not args.batch_root:
        args.batch_root = str(infer_batch_root(args.run_id))
        print("[infer] batch-root = " + args.batch_root)
    if not args.batch_id:
        try:
            args.batch_id = batch_id_from_run_id(args.run_id)
            print("[infer] batch-id = " + args.batch_id)
        except IntegrationError:
            args.batch_id = batch_id_from_root(args.batch_root)
    if not args.machine_id:
        args.machine_id = machine_from_run(args.run_id)
        print("[infer] machine-id = " + args.machine_id)

    batches_path = Path(args.batches)
    if not batches_path.exists():
        print("[init] batches.tsv 不存在，先从 YY2_info.txt 初始化")
        batches = init_batches_from_source(Path(args.source_info))
    else:
        batches = load_batches(batches_path)

    new_row = {
        "batch_id": args.batch_id or batch_id_from_root(args.batch_root),
        "runID": args.run_id,
        "machineID": args.machine_id or machine_from_run(args.run_id),
        "version": args.version,
        "batch_root": args.batch_root,
        "projects": ",".join(args.projects),
        "status": "done",
    }
    pending_batches = merge_batch(batches, new_row)
    issues = check_batch(args, pending_batches)
    run_work = Path(args.reuse_work) if args.reuse_work else WORK_DIR / ("dry-run-" + now_tag() if args.dry_run else "run-" + now_tag())
    run_work.mkdir(parents=True, exist_ok=True)

    if issues:
        print("[check] 发现以下问题:")
        for issue in issues:
            print("  - " + issue)
        if args.strict:
            raise IntegrationError("严格模式下检查未通过")

    if args.dry_run:
        print("[dry-run] 不会写入 batches.tsv 或 summary 正式目录")
    else:
        save_batches(batches_path, pending_batches)
        print("[batches] 已更新 " + str(batches_path))

    if not args.skip_stats:
        run_existing_project_scripts(args, active_batches(pending_batches, args.projects), run_work)

    if not args.dry_run:
        args.write = True
    counts = render_project_outputs(args, active_batches(pending_batches, args.projects), run_work)
    report_counts(counts, run_work)


def command_rebuild(args):
    ensure_dirs()
    batches = load_batches(Path(args.batches))
    if not batches:
        raise IntegrationError("batches.tsv 为空，请先 init-batches 或 add-batch")
    run_work = Path(args.reuse_work) if args.reuse_work else WORK_DIR / ("dry-run-" + now_tag() if args.dry_run else "run-" + now_tag())
    run_work.mkdir(parents=True, exist_ok=True)
    if not args.skip_stats:
        run_existing_project_scripts(args, active_batches(batches, args.projects), run_work)
    if not args.dry_run:
        args.write = True
    counts = render_project_outputs(args, active_batches(batches, args.projects), run_work)
    report_counts(counts, run_work)


def report_counts(counts, run_work):
    print("[summary] 输出预览/日志目录: " + str(run_work))
    for (project, key), count in counts.items():
        print("  %s %s: %s rows" % (project, key, count))


def build_parser():
    parser = argparse.ArgumentParser(description="Integrate YY2 PE batch summaries.")
    sub = parser.add_subparsers(dest="command", required=True)

    common_paths = argparse.ArgumentParser(add_help=False)
    common_paths.add_argument("--batches", default=str(CONFIG_DIR / "batches.tsv"))
    common_paths.add_argument("--source-info", default=str(DEFAULT_SOURCE_INFO))
    common_paths.add_argument("--summary-dir", default=str(DEFAULT_SUMMARY_DIR))
    common_paths.add_argument("--history-summary-dir", default=str(DEFAULT_SUMMARY_DIR))
    common_paths.add_argument("--pa-pos", default=str(DEFAULT_KB_FILE))
    common_paths.add_argument("--pb-pos", default=str(DEFAULT_PB_POS))
    common_paths.add_argument("--pd-pos", default=None)
    common_paths.add_argument("--dry-run", action="store_true")
    common_paths.add_argument("--write", action="store_true")
    common_paths.add_argument("--skip-stats", action="store_true", help="reuse existing files in work dir when iterating render code")
    common_paths.add_argument("--reuse-work", default="", help="existing integrate/work run directory to render from")
    common_paths.add_argument("--preview-files", action="store_true", help="write preview CSV files during dry-run")
    common_paths.add_argument("--strict", action="store_true")
    common_paths.add_argument(
        "--history-mode",
        choices=["summary", "raw"],
        default="summary",
        help="summary: reuse existing summary CSVs as historical data and compute only target batch; raw: recompute from all raw result directories",
    )
    common_paths.add_argument(
        "--existing-target-source",
        choices=["computed", "history", "fallback"],
        default="computed",
        help=(
            "computed: replace target rows with newly computed rows; "
            "history: if the target run already exists in history summary CSVs, keep those rows; "
            "fallback: use history target rows only when computed target rows are absent"
        ),
    )
    common_paths.add_argument(
        "--dp-scope",
        choices=["target", "all", "skip"],
        default="target",
        help="target: only compute depth table for --run-id and merge into existing CSV; all: full-history depth recompute; skip: do not build depth table",
    )

    p_init = sub.add_parser("init-batches", parents=[common_paths])
    p_init.set_defaults(func=command_init_batches)

    p_add = sub.add_parser("add-batch", parents=[common_paths])
    p_add.add_argument("--batch-root", default="", help="optional; inferred from --run-id SKII number when omitted")
    p_add.add_argument("--batch-id")
    p_add.add_argument("--run-id", required=True)
    p_add.add_argument("--machine-id")
    p_add.add_argument("--version", default="V26")
    p_add.add_argument("--projects", type=parse_projects, default=list(DEFAULT_PROJECTS))
    p_add.set_defaults(func=command_add_batch)

    p_rebuild = sub.add_parser("rebuild", parents=[common_paths])
    p_rebuild.add_argument("--batch-root", default="")
    p_rebuild.add_argument("--run-id", default="")
    p_rebuild.add_argument("--machine-id", default="")
    p_rebuild.add_argument("--version", default="")
    p_rebuild.add_argument("--projects", type=parse_projects, default=list(DEFAULT_PROJECTS))
    p_rebuild.set_defaults(func=command_rebuild)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except subprocess.CalledProcessError as exc:
        print("[error] 命令执行失败: %s" % exc, file=sys.stderr)
        return exc.returncode or 1
    except IntegrationError as exc:
        print("[error] %s" % exc, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
