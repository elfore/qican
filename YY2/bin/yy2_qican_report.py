#!/usr/bin/env python3
"""
Build a YY2 qican PA/PB/PD batch QC comparison report.
python3 yy2_qican_report.py \
  --batchid 260922153603_B183_SKII-JBJC-YY2-260922153621 \
[  --runid SKII \]
  --outdir yy2_report_out
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
from pathlib import Path
from statistics import mean
from typing import Any


DEFAULT_BASE = Path(
    "/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/"
    "01.BK/02.Analysis/04.new_qican/S100/PE150/summary/YY2"
)

PROJECTS = {
    "PA": {
        "file": "02_质控汇总.csv",
        "run_columns": ["runid"],
        "metrics": ["raw_reads", "clean_reads", "effective_rate"],
        "read_metrics": {"raw_reads", "clean_reads"},
        "qc_column": "质控",
    },
    "PB": {
        "file": "02_质控汇总.csv",
        "run_columns": ["runinfo"],
        "metrics": ["effective_rate", "on_target", "uniformity"],
        "read_metrics": set(),
    },
    "PD": {
        "file": "02_质控汇总.csv",
        "run_columns": ["runID"],
        "metrics": ["effective_rate", "on_target", "uniformity"],
        "read_metrics": set(),
    },
}

SAMPLE_COLUMNS = ["sample", "sample_id", "sampleID", "样本", "样本名", "样本编号"]


def read_csv(path: Path) -> list[dict[str, str]]:
    for encoding in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            with path.open("r", encoding=encoding, newline="") as handle:
                return list(csv.DictReader(handle))
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("csv", b"", 0, 1, f"Could not decode {path}")


def normalize_headers(row: dict[str, str]) -> dict[str, str]:
    return {(key or "").strip(): value for key, value in row.items()}


def first_existing(headers: list[str], candidates: list[str]) -> str | None:
    lowered = {header.lower(): header for header in headers}
    for candidate in candidates:
        if candidate in headers:
            return candidate
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def parse_number(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if not text or text in {"-", "NA", "N/A", "nan", "NaN"}:
        return None
    multiplier = 1.0
    if text.endswith("%"):
        text = text[:-1]
    if text.lower().endswith("m"):
        multiplier = 1_000_000.0
        text = text[:-1]
    elif text.lower().endswith("k"):
        multiplier = 1_000.0
        text = text[:-1]
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    if not match:
        return None
    return float(match.group(0)) * multiplier


def infer_batch_key(run_value: str, runid: str | None, batchid: str | None) -> str:
    if batchid:
        return batchid
    if not runid:
        return run_value
    if runid in run_value:
        start = run_value.find(runid)
        prefix = run_value[:start].rstrip("_- ")
        suffix = run_value[start:]
        parts = suffix.split()
        suffix = parts[0] if parts else suffix
        if prefix:
            return f"{prefix}_{suffix}"
        return suffix
    return run_value


def format_metric(metric: str, value: float | None, read_metrics: set[str]) -> str:
    if value is None:
        return "NA"
    if metric in read_metrics:
        return f"{value / 1_000_000:.2f}M"
    return f"{value:.2f}%"


def format_diff(metric: str, diff: float | None, read_metrics: set[str]) -> str:
    if diff is None:
        return "NA"
    sign = "+" if diff > 0 else ""
    if metric in read_metrics:
        return f"{sign}{diff / 1_000_000:.2f}M"
    return f"{sign}{diff:.2f}个百分点"


def trend(diff: float | None) -> str:
    if diff is None:
        return "NA"
    if abs(diff) < 1e-9:
        return "持平"
    return "上升" if diff > 0 else "下降"


def mean_for(rows: list[dict[str, str]], metric: str, *, percentage: bool) -> float | None:
    values = [parse_number(row.get(metric)) for row in rows]
    numeric = [value for value in values if value is not None]
    if percentage and numeric and max(abs(value) for value in numeric) <= 1:
        numeric = [value * 100 for value in numeric]
    return mean(numeric) if numeric else None


def row_matches(run_value: str, runid: str | None, batchid: str | None) -> bool:
    if batchid:
        return batchid in run_value
    if runid:
        return runid in run_value
    return False


def failed_pa_samples(rows: list[dict[str, str]], qc_column: str, run_column: str) -> list[str]:
    failed: list[str] = []
    for idx, row in enumerate(rows, start=1):
        qc = str(row.get(qc_column, "")).strip()
        if qc and qc == "合格":
            continue
        if not qc:
            continue
        sample_column = first_existing(list(row.keys()), SAMPLE_COLUMNS)
        sample = row.get(sample_column, "") if sample_column else ""
        if not sample:
            sample = f"第{idx}行({row.get(run_column, '')})"
        failed.append(f"{sample}: {qc}")
    return failed


def analyze_project(
    base: Path,
    project: str,
    config: dict[str, Any],
    runid: str | None,
    batchid: str | None,
) -> dict[str, Any]:
    path = base / project / config["file"]
    rows = [normalize_headers(row) for row in read_csv(path)]
    if not rows:
        raise ValueError(f"{path} is empty")

    headers = list(rows[0].keys())
    run_column = first_existing(headers, config["run_columns"]) or headers[0]
    matched = [row for row in rows if row_matches(str(row.get(run_column, "")), runid, batchid)]
    if not matched:
        query = batchid or runid or ""
        raise ValueError(f"No rows matched {query!r} in {path} column {run_column}")

    inferred_keys = sorted(
        {infer_batch_key(str(row.get(run_column, "")), runid, batchid) for row in matched}
    )
    if not batchid and runid and len(inferred_keys) > 1:
        return {
            "project": project,
            "ambiguous": True,
            "candidate_batches": inferred_keys,
            "run_column": run_column,
        }

    matched_ids = {id(row) for row in matched}
    history = [row for row in rows if id(row) not in matched_ids]
    read_metrics = set(config.get("read_metrics", set()))
    metrics = []
    for metric in config["metrics"]:
        is_percentage = metric not in read_metrics
        current_mean = mean_for(matched, metric, percentage=is_percentage)
        historical_mean = mean_for(history, metric, percentage=is_percentage)
        diff = None
        if current_mean is not None and historical_mean is not None:
            diff = current_mean - historical_mean
        metrics.append(
            {
                "metric": metric,
                "current_mean": current_mean,
                "historical_mean": historical_mean,
                "difference": diff,
                "current_display": format_metric(metric, current_mean, read_metrics),
                "historical_display": format_metric(metric, historical_mean, read_metrics),
                "difference_display": format_diff(metric, diff, read_metrics),
                "trend": trend(diff),
            }
        )

    result: dict[str, Any] = {
        "project": project,
        "path": str(path),
        "run_column": run_column,
        "matched_rows": len(matched),
        "history_rows": len(history),
        "batch_key": inferred_keys[0],
        "metrics": metrics,
    }

    if project == "PA" and config.get("qc_column") in headers:
        failed = failed_pa_samples(matched, config["qc_column"], run_column)
        result["pa_qc"] = {
            "all_passed": not failed,
            "failed_samples": failed,
            "message": "样本质控指标合格。" if not failed else "存在样本质控指标不合格。",
        }

    return result


def render_markdown(report: dict[str, Any]) -> str:
    lines = [
        f"# YY2测序评估批次质控报告 - {report['batch_label']}",
        "",
        "## 一、批次信息",
        "",
        f"Batch ID：{report['batch_label']}",
        "",
    ]
    section_names = {"PA": "二、PA项目", "PB": "三、PB项目", "PD": "四、PD项目"}
    for project in ("PA", "PB", "PD"):
        item = report["projects"][project]
        lines.extend(
            [
                f"## {section_names[project]}",
                "",
                "| 指标 | 本批次均值 | 历史批次均值 | 差值 | 趋势 |",
                "| --- | --- | --- | --- | --- |",
            ]
        )
        for metric in item["metrics"]:
            lines.append(
                "| {metric} | {current_display} | {historical_display} | "
                "{difference_display} | {trend} |".format(**metric)
            )
        lines.append("")
        if project == "PA" and "pa_qc" in item:
            if item["pa_qc"]["all_passed"]:
                lines.append("样本质控指标合格。")
            else:
                lines.append("质控不合格样本：" + "；".join(item["pa_qc"]["failed_samples"]))
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_html(report: dict[str, Any]) -> str:
    parts = [
        f"<h1>YY2测序评估批次质控报告 - {html.escape(report['batch_label'])}</h1>",
        "<h2>一、批次信息</h2>",
        f"<p><strong>Batch ID：</strong>{html.escape(report['batch_label'])}</p>",
    ]
    section_names = {"PA": "二、PA项目", "PB": "三、PB项目", "PD": "四、PD项目"}
    for project in ("PA", "PB", "PD"):
        item = report["projects"][project]
        parts.append(f"<h2>{section_names[project]}</h2>")
        parts.append("<table><thead><tr>")
        for header in ["指标", "本批次均值", "历史批次均值", "差值", "趋势"]:
            parts.append(f"<th>{header}</th>")
        parts.append("</tr></thead><tbody>")
        for metric in item["metrics"]:
            parts.append("<tr>")
            for key in ["metric", "current_display", "historical_display", "difference_display", "trend"]:
                parts.append(f"<td>{html.escape(str(metric[key]))}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table>")
        if project == "PA" and "pa_qc" in item:
            if item["pa_qc"]["all_passed"]:
                parts.append("<p>样本质控指标合格。</p>")
            else:
                failed = "；".join(html.escape(x) for x in item["pa_qc"]["failed_samples"])
                parts.append(f"<p>质控不合格样本：{failed}</p>")
    return "\n".join(parts) + "\n"


def build_report(base: Path, runid: str | None, batchid: str | None) -> dict[str, Any]:
    if not runid and not batchid:
        raise ValueError("Provide --runid or --batchid")
    projects = {
        project: analyze_project(base, project, config, runid, batchid)
        for project, config in PROJECTS.items()
    }
    ambiguous = {name: data for name, data in projects.items() if data.get("ambiguous")}
    if ambiguous:
        candidates = sorted(
            {
                candidate
                for data in ambiguous.values()
                for candidate in data.get("candidate_batches", [])
            }
        )
        raise SystemExit(
            "Run id is ambiguous. Re-run with --batchid. Candidates:\n"
            + "\n".join(f"- {candidate}" for candidate in candidates)
        )
    batch_label = batchid or next(iter(projects.values()))["batch_key"]
    return {"batch_label": batch_label, "query": {"runid": runid, "batchid": batchid}, "projects": projects}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE, help="YY2 summary directory")
    parser.add_argument("--runid", help="Short run id, for example SKII")
    parser.add_argument("--batchid", help="Full batch id")
    parser.add_argument("--outdir", type=Path, default=Path("yy2_qican_report_out"))
    args = parser.parse_args()

    report = build_report(args.base, args.runid, args.batchid)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (args.outdir / "report.md").write_text(render_markdown(report), encoding="utf-8")
    (args.outdir / "report.html").write_text(render_html(report), encoding="utf-8")
    print(args.outdir.resolve())


if __name__ == "__main__":
    main()
