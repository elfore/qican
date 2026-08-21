#!/usr/bin/env python3
"""
Workflow wrapper for YY2 qican PE150 PA/PB/PD summary integration.

This script orchestrates qican_integrate.py without duplicating its
calculation logic. It adds runID inference, prechecks, run modes, validation,
and a concise report.
"""

import argparse
import csv
import datetime as _dt
import shutil
import subprocess
import sys
from pathlib import Path

import qican_integrate as qi


SCRIPT_DIR = Path(__file__).resolve().parent
INTEGRATOR = SCRIPT_DIR / "qican_integrate.py"
BASE_DIR = SCRIPT_DIR.parent
REPORT_DIR = BASE_DIR / "reports"
LOG_DIR = BASE_DIR / "logs"


class WorkflowError(RuntimeError):
    pass


def now_tag():
    return _dt.datetime.now().strftime("%Y%m%d-%H%M%S")


def parse_projects(value):
    return qi.parse_projects(value)


def read_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.reader(handle))


def max_cols(rows):
    return max((len(row) for row in rows), default=0)


def empty_columns(rows):
    cols = max_cols(rows)
    return [
        idx + 1
        for idx in range(cols)
        if all(not (idx < len(row) and str(row[idx]).strip()) for row in rows)
    ]


def target_rows(path, run_id):
    rows = read_csv(path)
    if not rows:
        return []
    return [row for row in rows[1:] if row and row[0] == run_id]


def matrix_has_target(path, run_id):
    for row in read_csv(path):
        if len(row) >= 3 and row[1] == "一致率":
            if any(value.startswith(run_id + "_") for value in row[2:]):
                return True
    return False


def infer_run(run_id):
    batch_id = qi.batch_id_from_run_id(run_id)
    batch_root = qi.infer_batch_root(run_id)
    machine_id = qi.machine_from_run(run_id)
    return {
        "run_id": run_id,
        "batch_id": batch_id,
        "batch_root": str(batch_root),
        "machine_id": machine_id,
    }


def precheck(run_info, projects):
    batch_root = Path(run_info["batch_root"])
    issues = []
    warnings = []
    if not batch_root.exists():
        issues.append("batch_root missing: " + str(batch_root))
    for project in projects:
        project_dir = batch_root / project
        if not project_dir.exists():
            issues.append("%s dir missing: %s" % (project, project_dir))
    if "PB" in projects and (batch_root / "PB").exists() and not (batch_root / "PB" / "result").exists():
        warnings.append("PB/result missing; PB summary may be incomplete")
    if "PD" in projects and not qi.find_pd_pos(batch_root):
        warnings.append("PD result/mutation/*final_mut.txt not found; PD 01/03/04 may be incomplete")
    return issues, warnings


def make_temp_root(mode, batch_id):
    root = Path("/tmp") / ("qican_workflow_%s_%s_%s" % (mode, batch_id, now_tag()))
    root.mkdir(parents=True, exist_ok=True)
    return root


def prepare_temp_batches(temp_root, batches_path):
    dst = temp_root / "batches.tsv"
    src = Path(batches_path)
    if src.exists():
        shutil.copy2(str(src), str(dst))
    return dst


def stream_command(cmd, log_handle):
    log_handle.write("[cmd] " + " ".join(map(str, cmd)) + "\n")
    log_handle.flush()
    proc = subprocess.Popen(
        [str(x) for x in cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line, end="")
        log_handle.write(line)
    return proc.wait()


def build_integrate_cmd(args, run_info, temp_root):
    cmd = [
        sys.executable,
        INTEGRATOR,
        "add-batch",
        "--run-id",
        run_info["run_id"],
        "--history-mode",
        "summary",
        "--dp-scope",
        args.dp_scope,
        "--version",
        args.version,
        "--projects",
        ",".join(args.projects),
        "--history-summary-dir",
        args.history_summary_dir,
    ]
    if args.mode == "dry-run":
        cmd.append("--dry-run")
        if temp_root:
            cmd.extend(["--reuse-work", str(temp_root / "work")])
    elif args.mode == "regression":
        if not temp_root:
            raise WorkflowError("regression mode requires a temporary root")
        temp_batches = prepare_temp_batches(temp_root, args.batches)
        cmd.extend(
            [
                "--batches",
                str(temp_batches),
                "--summary-dir",
                str(temp_root / "summary"),
                "--reuse-work",
                str(temp_root / "work"),
                "--existing-target-source",
                "history",
            ]
        )
    else:
        cmd.extend(["--summary-dir", args.summary_dir])
    if args.skip_stats:
        cmd.append("--skip-stats")
    if args.preview_files:
        cmd.append("--preview-files")
    if args.strict:
        cmd.append("--strict")
    return cmd


def validation_summary(summary_dir, run_ids, projects):
    summary = Path(summary_dir)
    result = {
        "summary_dir": str(summary),
        "empty_col_files": [],
        "missing_target_rows": [],
        "missing_target_matrix_columns": [],
    }
    for project in projects:
        project_dir = summary / project
        if not project_dir.exists():
            result["missing_target_rows"].append((project, "PROJECT_DIR_MISSING"))
            continue
        for csv_path in sorted(project_dir.glob("*.csv")):
            rows = read_csv(csv_path)
            empty = empty_columns(rows)
            if empty:
                result["empty_col_files"].append((project, csv_path.name, empty[0], empty[-1], len(empty)))
    for run_id in run_ids:
        for project in projects:
            for key in ("01", "02", "03", "04"):
                filename = qi.SUMMARY_FILENAMES.get((project, key))
                if not filename:
                    continue
                path = summary / project / filename
                if path.exists() and not target_rows(path, run_id):
                    result["missing_target_rows"].append((run_id, project, key, filename))
        for project in ("PB", "PD"):
            if project not in projects:
                continue
            for key in ("05", "06", "07"):
                filename = qi.SUMMARY_FILENAMES.get((project, key))
                path = summary / project / filename
                if path.exists() and not matrix_has_target(path, run_id):
                    result["missing_target_matrix_columns"].append((run_id, project, key, filename))
    return result


def write_report(report_path, entries, validations):
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write("YY2 qican summary workflow report\n")
        handle.write("generated_at\t%s\n\n" % _dt.datetime.now().isoformat(timespec="seconds"))
        for entry in entries:
            handle.write("[%s]\n" % entry["run_id"])
            for key in ("mode", "batch_id", "machine_id", "batch_root", "exit_code", "log"):
                handle.write("%s\t%s\n" % (key, entry.get(key, "")))
            if entry.get("issues"):
                handle.write("issues\t%s\n" % "; ".join(entry["issues"]))
            if entry.get("warnings"):
                handle.write("warnings\t%s\n" % "; ".join(entry["warnings"]))
            handle.write("\n")
        for validation in validations:
            handle.write("[validation]\n")
            handle.write("summary_dir\t%s\n" % validation["summary_dir"])
            handle.write("empty_col_files\t%d\t%s\n" % (len(validation["empty_col_files"]), validation["empty_col_files"]))
            handle.write(
                "missing_target_rows\t%d\t%s\n"
                % (len(validation["missing_target_rows"]), validation["missing_target_rows"])
            )
            handle.write(
                "missing_target_matrix_columns\t%d\t%s\n"
                % (len(validation["missing_target_matrix_columns"]), validation["missing_target_matrix_columns"])
            )
            handle.write("\n")


def run_workflow(args):
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    tag = now_tag()
    entries = []
    validation_targets = []
    successful_write_run_ids = []

    for run_id in args.run_ids:
        run_info = infer_run(run_id)
        issues, warnings = precheck(run_info, args.projects)
        temp_root = make_temp_root(args.mode, run_info["batch_id"]) if args.mode in ("dry-run", "regression") else None
        log_path = LOG_DIR / ("workflow_%s_%s_%s.log" % (args.mode, run_info["batch_id"], tag))
        print("[workflow] %s -> %s %s" % (run_id, run_info["batch_id"], run_info["machine_id"]))
        print("[workflow] batch-root: " + run_info["batch_root"])
        for warning in warnings:
            print("[warn] " + warning)
        if issues:
            for issue in issues:
                print("[issue] " + issue)
            if args.strict:
                raise WorkflowError("precheck failed for " + run_id)

        cmd = build_integrate_cmd(args, run_info, temp_root)
        with open(log_path, "w", encoding="utf-8") as log_handle:
            log_handle.write("[precheck] issues=%s warnings=%s\n" % (issues, warnings))
            exit_code = stream_command(cmd, log_handle)

        entry = {
            "mode": args.mode,
            "run_id": run_id,
            "batch_id": run_info["batch_id"],
            "machine_id": run_info["machine_id"],
            "batch_root": run_info["batch_root"],
            "exit_code": exit_code,
            "issues": issues,
            "warnings": warnings,
            "log": str(log_path),
        }
        entries.append(entry)
        if exit_code != 0 and not args.keep_going:
            break
        if args.mode == "write":
            successful_write_run_ids.append(run_id)
            validation_targets.append((Path(args.summary_dir), tuple(successful_write_run_ids)))
        elif args.mode == "regression" and temp_root:
            validation_targets.append((temp_root / "summary", (run_id,)))

    validations = []
    if args.mode in ("write", "regression"):
        seen = set()
        for summary_dir, run_ids in validation_targets:
            key = (str(summary_dir), run_ids)
            if key in seen:
                continue
            seen.add(key)
            validations.append(validation_summary(summary_dir, run_ids, args.projects))

    report_path = REPORT_DIR / ("workflow_%s_%s.txt" % (args.mode, tag))
    write_report(report_path, entries, validations)
    print("[workflow] report: " + str(report_path))
    for validation in validations:
        print("[validation] summary_dir: " + validation["summary_dir"])
        print("[validation] empty_col_files: %d %s" % (len(validation["empty_col_files"]), validation["empty_col_files"]))
        print(
            "[validation] missing_target_rows: %d %s"
            % (len(validation["missing_target_rows"]), validation["missing_target_rows"])
        )
        print(
            "[validation] missing_target_matrix_columns: %d %s"
            % (len(validation["missing_target_matrix_columns"]), validation["missing_target_matrix_columns"])
        )
    return max((entry["exit_code"] for entry in entries), default=0)


def build_parser():
    parser = argparse.ArgumentParser(description="YY2 qican PE150 summary workflow.")
    parser.add_argument("run_ids", nargs="+", help="original batchid/runID strings")
    parser.add_argument("--mode", choices=["dry-run", "write", "regression"], default="dry-run")
    parser.add_argument("--version", default="V26")
    parser.add_argument("--projects", type=parse_projects, default=list(qi.DEFAULT_PROJECTS))
    parser.add_argument("--dp-scope", choices=["target", "all", "skip"], default="target")
    parser.add_argument("--batches", default=str(qi.CONFIG_DIR / "batches.tsv"))
    parser.add_argument("--summary-dir", default=str(qi.DEFAULT_SUMMARY_DIR))
    parser.add_argument("--history-summary-dir", default=str(qi.DEFAULT_SUMMARY_DIR))
    parser.add_argument("--skip-stats", action="store_true")
    parser.add_argument("--preview-files", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return run_workflow(args)
    except (WorkflowError, qi.IntegrationError) as exc:
        print("[error] " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
