#!/usr/bin/env python3
"""Export YY1 mNGS Feishu abundance history sheets to long TSV.

The script is intended to run on the cluster after sourcing the Codex
environment that provides lark-cli. It reads Feishu sheet display values via
`lark-cli sheets +csv-get`, so percentage cells are parsed from the visible
values instead of typed table values.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable


DEFAULT_URL = "https://cygnusbio.feishu.cn/wiki/YwevwkHrSiTUlJkPmuLc2COznWb?sheet=b6caGo"
DEFAULT_LARK_ENV = "/mnt/gpfs1/Users/yangjinxurong/software/codex/env.sh"
OUTPUT_COLUMNS = [
    "sample",
    "batch",
    "taxid",
    "species",
    "gc_pct",
    "expected_pct",
    "abundance_pct",
]
SKII15404_EXPECTED_COUNTS = {
    "CN002341-S01-D01-L02-UDI1": 10,
    "CN002386-S01-D01-L01-UDI2": 11,
    "CN002387-S01-D01-L01-UDI3": 10,
    "CN002388-S01-D01-L01-UDI4": 10,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export Feishu YY1 mNGS abundance history to long TSV."
    )
    parser.add_argument("--url", default=DEFAULT_URL, help="Feishu wiki/sheet URL.")
    parser.add_argument(
        "--spreadsheet-token",
        default=None,
        help="Use a spreadsheet token directly instead of resolving --url.",
    )
    parser.add_argument(
        "--out",
        default="feishu_abundance_history_long.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--current-batch",
        default=None,
        help="Optional batch name to validate and summarize, e.g. SKII15404.",
    )
    parser.add_argument(
        "--as",
        dest="identity",
        default="user",
        choices=("user", "bot"),
        help="lark-cli identity to use. Default: user.",
    )
    parser.add_argument(
        "--as-user",
        action="store_const",
        const="user",
        dest="identity",
        help="Explicitly use lark-cli user identity.",
    )
    parser.add_argument(
        "--lark-env",
        default=DEFAULT_LARK_ENV,
        help="Shell env file to source before running lark-cli. Use '' to disable.",
    )
    parser.add_argument(
        "--sample-prefix",
        default="CN",
        help="Only sheets whose names start with this prefix are exported.",
    )
    return parser.parse_args()


def run_lark(args: list[str], env_file: str | None, identity: str) -> dict[str, Any]:
    full_args = ["lark-cli", *args, "--as", identity, "--format", "json"]
    env = os.environ.copy()
    env["LARKSUITE_CLI_NO_UPDATE_NOTIFIER"] = "1"
    env["LARKSUITE_CLI_NO_SKILLS_NOTIFIER"] = "1"

    if env_file:
        quoted = " ".join(shlex.quote(part) for part in full_args)
        cmd: list[str] = [
            "bash",
            "-lc",
            f"source {shlex.quote(env_file)} >/dev/null 2>&1; {quoted}",
        ]
    else:
        cmd = full_args

    proc = subprocess.run(cmd, text=True, capture_output=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(
            "lark-cli failed with exit code "
            f"{proc.returncode}\nCOMMAND: {' '.join(shlex.quote(x) for x in full_args)}"
            f"\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )

    text = proc.stdout.strip()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "lark-cli did not return valid JSON\n"
            f"COMMAND: {' '.join(shlex.quote(x) for x in full_args)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        ) from exc

    if not payload.get("ok", False):
        raise RuntimeError(
            "lark-cli returned ok=false\n"
            f"COMMAND: {' '.join(shlex.quote(x) for x in full_args)}\n"
            f"RESPONSE:\n{json.dumps(payload, ensure_ascii=False, indent=2)}"
        )
    return payload


def resolve_spreadsheet_token(
    url: str, spreadsheet_token: str | None, env_file: str | None, identity: str
) -> tuple[list[str], str]:
    if spreadsheet_token:
        return ["--spreadsheet-token", spreadsheet_token], spreadsheet_token

    inspect = run_lark(["drive", "+inspect", "--url", url], env_file, identity)
    data = inspect.get("data", {})
    doc_type = data.get("type") or data.get("obj_type")
    token = data.get("token") or data.get("obj_token")
    if doc_type != "sheet" or not token:
        raise RuntimeError(
            f"Input URL did not resolve to a sheet. type={doc_type!r}, token={token!r}"
        )
    return ["--spreadsheet-token", str(token)], str(token)


def excel_col_name(index_1based: int) -> str:
    if index_1based < 1:
        raise ValueError(f"Column index must be >= 1: {index_1based}")
    chars: list[str] = []
    n = index_1based
    while n:
        n, remainder = divmod(n - 1, 26)
        chars.append(chr(ord("A") + remainder))
    return "".join(reversed(chars))


def strip_row_prefix(line: str) -> str:
    return re.sub(r"^\[row=\d+\]\s?", "", line)


def parse_annotated_csv(annotated_csv: str) -> list[list[str]]:
    cleaned = "\n".join(strip_row_prefix(line) for line in annotated_csv.splitlines())
    return list(csv.reader(io.StringIO(cleaned)))


def normalize_header(value: str) -> str:
    return re.sub(r"\s+", "", value.strip().lower())


def is_meta_header(header: str) -> bool:
    norm = normalize_header(header)
    return norm in {
        "species",
        "物种",
        "expeted_speices",
        "taxid",
        "tax_id",
        "taxonid",
        "gc",
        "gc%",
        "gc_pct",
        "预期",
        "预期比例",
        "expected",
        "expected%",
        "expected_pct",
        "maxdpl",
        "maxhpl",
    }


def find_column(headers: list[str], candidates: Iterable[str]) -> int | None:
    wanted = {normalize_header(item) for item in candidates}
    for idx, header in enumerate(headers):
        if normalize_header(header) in wanted:
            return idx
    return None


def cell(row: list[str], idx: int | None) -> str:
    if idx is None or idx >= len(row):
        return ""
    return row[idx].strip()


def parse_percent_or_number(value: str) -> str:
    text = value.strip()
    if not text:
        return ""
    text = text.replace(",", "")
    if text.endswith("%"):
        text = text[:-1].strip()
    try:
        number = Decimal(text)
    except InvalidOperation:
        return value.strip()
    return format(number.normalize(), "f")


def sheet_title(sheet: dict[str, Any]) -> str:
    return str(sheet.get("sheet_name") or sheet.get("title") or sheet.get("name") or "")


def sheet_id(sheet: dict[str, Any]) -> str:
    return str(sheet.get("sheet_id") or sheet.get("id") or "")


def sheet_dimension(sheet: dict[str, Any], key: str, default: int) -> int:
    value = sheet.get(key)
    if value is None:
        value = sheet.get("grid_properties", {}).get(key)
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def export_sheet(
    locator: list[str],
    sheet: dict[str, Any],
    env_file: str | None,
    identity: str,
) -> list[dict[str, str]]:
    sample = sheet_title(sheet)
    sid = sheet_id(sheet)
    row_count = sheet_dimension(sheet, "row_count", 1000)
    col_count = sheet_dimension(sheet, "column_count", 200)
    end_cell = f"{excel_col_name(col_count)}{row_count}"

    csv_payload = run_lark(
        [
            "sheets",
            "+csv-get",
            *locator,
            "--sheet-id",
            sid,
            "--range",
            f"A1:{end_cell}",
        ],
        env_file,
        identity,
    )
    data = csv_payload.get("data", {})
    annotated = data.get("annotated_csv") or data.get("csv") or ""
    if not annotated:
        raise RuntimeError(f"No annotated_csv returned for sheet {sample} ({sid})")

    rows = parse_annotated_csv(annotated)
    if not rows:
        return []

    headers = [item.strip() for item in rows[0]]
    species_idx = find_column(headers, ["species", "物种"])
    if species_idx is None:
        species_idx = find_column(headers, ["expeted_speices"])
    if species_idx is None:
        raise RuntimeError(f"Sheet {sample} has no species column")
    taxid_idx = find_column(headers, ["taxid", "tax_id", "taxonid"])
    gc_idx = find_column(headers, ["GC", "GC%", "gc_pct"])
    expected_idx = find_column(
        headers, ["预期", "预期比例", "expected", "expected%", "expected_pct"]
    )

    batch_columns = [
        (idx, header.strip())
        for idx, header in enumerate(headers)
        if header.strip() and not is_meta_header(header)
    ]

    out_rows: list[dict[str, str]] = []
    for row in rows[1:]:
        species = cell(row, species_idx)
        if not species:
            continue
        taxid = cell(row, taxid_idx) or species
        gc_pct = parse_percent_or_number(cell(row, gc_idx))
        expected_pct = parse_percent_or_number(cell(row, expected_idx))
        for idx, batch in batch_columns:
            abundance = parse_percent_or_number(cell(row, idx))
            if not abundance:
                continue
            out_rows.append(
                {
                    "sample": sample,
                    "batch": batch,
                    "taxid": taxid,
                    "species": species,
                    "gc_pct": gc_pct,
                    "expected_pct": expected_pct,
                    "abundance_pct": abundance,
                }
            )
    return out_rows


def validate(rows: list[dict[str, str]], current_batch: str | None) -> None:
    samples = sorted({row["sample"] for row in rows})
    print(f"Exported {len(rows)} rows from {len(samples)} sample sheets.", file=sys.stderr)
    for sample in samples:
        print(f"  - {sample}", file=sys.stderr)

    if not current_batch:
        return

    current_rows = [row for row in rows if row["batch"] == current_batch]
    counts = Counter(row["sample"] for row in current_rows)
    print(f"{current_batch} row counts:", file=sys.stderr)
    for sample in samples:
        print(f"  - {sample}: {counts.get(sample, 0)}", file=sys.stderr)

    mismatches = {
        sample: (expected, counts.get(sample, 0))
        for sample, expected in SKII15404_EXPECTED_COUNTS.items()
        if counts.get(sample, 0) != expected
    }
    if mismatches:
        details = ", ".join(
            f"{sample} expected {exp}, got {got}"
            for sample, (exp, got) in mismatches.items()
        )
        raise RuntimeError(f"SKII15404 validation failed: {details}")

    key_rows = [
        row
        for row in current_rows
        if row["sample"] == "CN002387-S01-D01-L01-UDI3"
        and row["species"] == "Streptococcus pyogenes"
    ]
    if key_rows:
        print(
            "Check CN002387-S01-D01-L01-UDI3 / SKII15404 / "
            "Streptococcus pyogenes abundance_pct = "
            f"{key_rows[0]['abundance_pct']}",
            file=sys.stderr,
            )


def main() -> int:
    args = parse_args()
    env_file = args.lark_env or None

    locator, token = resolve_spreadsheet_token(
        args.url, args.spreadsheet_token, env_file, args.identity
    )
    print(f"Using spreadsheet token: {token}", file=sys.stderr)

    workbook = run_lark(["sheets", "+workbook-info", *locator], env_file, args.identity)
    sheets = workbook.get("data", {}).get("sheets", [])
    sample_sheets = [
        sheet for sheet in sheets if sheet_title(sheet).startswith(args.sample_prefix)
    ]
    if not sample_sheets:
        raise RuntimeError(f"No sheets start with prefix {args.sample_prefix!r}")

    rows: list[dict[str, str]] = []
    for sheet in sample_sheets:
        title = sheet_title(sheet)
        print(f"Reading {title} ...", file=sys.stderr)
        rows.extend(export_sheet(locator, sheet, env_file, args.identity))

    out_path = Path(args.out)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    validate(rows, args.current_batch)
    print(f"Wrote {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
