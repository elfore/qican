#!/usr/bin/env python3
"""Safely upsert one YY1 batch into the four mNGS abundance Feishu sheets.

Only the named batch column is changed. Existing values, formulas, layout and
formats are preserved. A missing species or an unknown library is a hard
failure before any write is attempted.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
import sys
from pathlib import Path


YY1_DIR = Path("/mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1")
LARK = "/mnt/gpfs1/Users/yangjinxurong/software/lark-cli/bin/lark-cli"
TOKEN = "MjihsiLrnh6UU5tc4cSc4PRSntc"
EXPECTED_LIBRARIES = {
    "CN002341-S01-D01-L02-UDI1",
    "CN002386-S01-D01-L01-UDI2",
    "CN002387-S01-D01-L01-UDI3",
    "CN002388-S01-D01-L01-UDI4",
}


class SyncError(RuntimeError):
    pass


def column_name(index: int) -> str:
    result = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        result = chr(65 + remainder) + result
    return result


def run_cli(args: list[str], stdin: str | None = None) -> dict:
    proc = subprocess.run([LARK, *args, "--format", "json"], input=stdin, text=True, capture_output=True)
    if proc.returncode:
        raise SyncError(f"lark-cli failed: {' '.join(args)}\n{proc.stderr}")
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise SyncError(f"Invalid lark-cli JSON: {proc.stdout}\n{proc.stderr}") from exc
    if not payload.get("ok"):
        raise SyncError(f"lark-cli returned ok=false: {json.dumps(payload, ensure_ascii=False)}")
    return payload


def manifest_row(batch: str) -> dict[str, str]:
    info = YY1_DIR / "bin/YY1_info.txt"
    with info.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[0].casefold() == batch.casefold():
                return {"batch": fields[0], "run_id": fields[1], "dir": fields[3]}
    raise SyncError(f"Batch {batch!r} is not registered in {info}")


def parse_source(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 2:
        raise SyncError(f"No abundance data in {path}")
    header = [value.casefold().strip() for value in rows[0]]
    has_taxid = "taxid" in header
    values: dict[str, str] = {}
    for row in rows[1:]:
        if len(row) < 2 or not row[0].strip():
            continue
        value = row[-1].strip() if has_taxid else row[1].strip()
        if not value:
            raise SyncError(f"Empty abundance for {row[0]!r} in {path}")
        if row[0] in values:
            raise SyncError(f"Duplicate species {row[0]!r} in {path}")
        values[row[0]] = value
    if not values:
        raise SyncError(f"No valid abundance rows in {path}")
    return values


def parse_annotated_csv(value: str) -> list[list[str]]:
    # csv.reader handles quoted commas; row prefixes are CLI metadata only.
    cleaned = "\n".join(re.sub(r"^\[row=\d+\]\s?", "", line) for line in value.splitlines())
    return list(csv.reader(io.StringIO(cleaned)))


def read_sheet(sheet: dict) -> list[list[str]]:
    end = f"{column_name(int(sheet['column_count']))}{int(sheet['row_count'])}"
    payload = run_cli([
        "sheets", "+csv-get", "--spreadsheet-token", TOKEN,
        "--sheet-id", sheet["sheet_id"], "--range", f"A1:{end}",
    ])
    rows = parse_annotated_csv(payload["data"].get("annotated_csv", ""))
    if len(rows) < 2 or not rows[0]:
        raise SyncError(f"Sheet {sheet['sheet_name']} has no usable header/data")
    return rows


def plan_one(batch: str, sheet: dict, source: dict[str, str]) -> tuple[int, list[tuple[int, str]], list[list[str]]]:
    rows = read_sheet(sheet)
    header = rows[0]
    try:
        col_index = header.index(batch)
    except ValueError:
        # The caller may create this trailing column only after every sheet has
        # passed preflight. Values will inherit the style of the prior column.
        col_index = len(header)
    species_rows: dict[str, int] = {}
    for row_index, row in enumerate(rows[1:], start=2):
        if row and row[0].strip():
            species_rows[row[0].strip()] = row_index
    missing = sorted(set(source) - set(species_rows))
    if missing:
        raise SyncError(
            f"{sheet['sheet_name']}: source species absent from Feishu template: {', '.join(missing)}"
        )
    updates = sorted((species_rows[species], value) for species, value in source.items())
    target_rows = [row for row, _ in updates]
    if target_rows != list(range(target_rows[0], target_rows[-1] + 1)):
        raise SyncError(f"{sheet['sheet_name']}: source species are not a contiguous template block")
    return col_index, updates, rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", help="e.g. BGISEQ_0707")
    parser.add_argument("--apply", action="store_true", help="Perform writes after full preflight; default is dry-run.")
    args = parser.parse_args()

    item = manifest_row(args.batch)
    abundance_dir = Path(item["dir"]) / "mNGS"
    source_files = {path.stem: parse_source(path) for path in abundance_dir.glob("CN*.txt")}
    missing = EXPECTED_LIBRARIES - set(source_files)
    unexpected = set(source_files) - EXPECTED_LIBRARIES
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing=" + ",".join(sorted(missing)))
        if unexpected:
            details.append("unknown=" + ",".join(sorted(unexpected)))
        raise SyncError("Unexpected mNGS library set: " + "; ".join(details))

    workbook = run_cli(["sheets", "+workbook-info", "--spreadsheet-token", TOKEN])
    sheets = {sheet["sheet_name"]: sheet for sheet in workbook["data"]["sheets"]}
    missing_sheets = EXPECTED_LIBRARIES - set(sheets)
    if missing_sheets:
        raise SyncError("Missing mNGS Feishu sheets: " + ", ".join(sorted(missing_sheets)))

    plans = []
    for library in sorted(EXPECTED_LIBRARIES):
        col_index, updates, rows = plan_one(item["batch"], sheets[library], source_files[library])
        plans.append((library, sheets[library], col_index, updates, rows))
        action = "update" if col_index < len(rows[0]) else "append"
        print(f"PLAN {library}: {action} column {column_name(col_index + 1)}; {len(source_files[library])} species")

    if not args.apply:
        print("DRY-RUN OK: no Feishu data changed")
        return 0

    for library, sheet, col_index, updates, rows in plans:
        target = column_name(col_index + 1)
        if col_index == len(rows[0]):
            # csv-put safely auto-expands at the sheet boundary.  The structure
            # endpoint rejects a position one column past the current boundary.
            run_cli([
                "sheets", "+csv-put", "--spreadsheet-token", TOKEN,
                "--sheet-id", sheet["sheet_id"], "--start-cell", f"{target}1", "--csv", "-",
            ], stdin=item["batch"] + "\n")
        csv_buffer = io.StringIO()
        csv.writer(csv_buffer, lineterminator="\n").writerows([[value] for _, value in updates])
        payload = csv_buffer.getvalue()
        run_cli([
            "sheets", "+csv-put", "--spreadsheet-token", TOKEN,
            "--sheet-id", sheet["sheet_id"], "--start-cell", f"{target}{updates[0][0]}", "--csv", "-",
        ], stdin=payload)
        # A newly appended column changes the sheet dimensions; re-read the
        # workbook metadata before validating the write-back.
        refreshed = run_cli(["sheets", "+workbook-info", "--spreadsheet-token", TOKEN])
        refreshed_sheet = next(item for item in refreshed["data"]["sheets"] if item["sheet_id"] == sheet["sheet_id"])
        verified = read_sheet(refreshed_sheet)
        actual = [verified[row - 1][col_index] if col_index < len(verified[row - 1]) else "" for row, _ in updates]
        if actual != [value for _, value in updates]:
            raise SyncError(f"Read-back validation failed for {library} / {target}")
        print(f"SYNCED {library}: {target} ({len(updates)} species)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SyncError as exc:
        print(f"SYNC FAILED: {exc}", file=sys.stderr)
        raise SystemExit(2)
