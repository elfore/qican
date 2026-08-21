#!/usr/bin/env python3
"""Build validated, versioned YY1 master CSV tables from assessment batches.

The files written by this program are the local source of truth.  A Feishu
syncer must consume these files only after this program exits successfully.
No Feishu write is performed here.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path("/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY1")
DEFAULT_MANIFEST = PROJECT_ROOT / "bin/YY1_info.txt"
DEFAULT_OUT = Path(
    "/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/"
    "02.Analysis/04.new_qican/S100/SE75/summary/YY1/csv"
)


class ValidationError(RuntimeError):
    """Raised before any master table is replaced."""


def read_tsv(path: Path, *, header: bool = False) -> tuple[list[str], list[list[str]]]:
    if not path.is_file():
        raise ValidationError(f"Missing required file: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if any(cell.strip() for cell in row)]
    if not rows:
        raise ValidationError(f"Empty required file: {path}")
    return (rows.pop(0), rows) if header else ([], rows)


def write_csv_atomic(path: Path, rows: list[list[str]]) -> dict[str, int]:
    if not rows:
        raise ValidationError(f"Refusing to replace {path}: no rows")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, delete=False) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows(rows)
        temp_path = Path(handle.name)
    os.replace(temp_path, path)
    return {"rows": len(rows), "columns": max(len(row) for row in rows)}


def parse_manifest(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    _, raw = read_tsv(path)
    for number, row in enumerate(raw, start=1):
        if len(row) < 4:
            raise ValidationError(f"{path}:{number}: expected 4 tab-separated fields")
        batch, run_id, read_type, directory = row[:4]
        if batch in seen:
            raise ValidationError(f"Duplicate batch in {path}: {batch}")
        directory_path = Path(directory)
        if not directory_path.is_dir():
            raise ValidationError(f"Batch {batch}: result directory does not exist: {directory}")
        seen.add(batch)
        rows.append({"batch": batch, "run_id": run_id, "read_type": read_type, "dir": directory})
    return rows


def existing_global_table(path: Path) -> list[list[str]]:
    _, rows = read_tsv(path)
    return rows


def result_file(root: Path, *relative_paths: str) -> Path:
    """Return the first existing layout-compatible result path."""
    candidates = [root / relative for relative in relative_paths]
    return next((path for path in candidates if path.exists()), candidates[0])


def result_dir(root: Path, *relative_paths: str) -> Path:
    """Return the first existing layout-compatible result directory."""
    candidates = [root / relative for relative in relative_paths]
    return next((path for path in candidates if path.is_dir()), candidates[0])


def lambda_summary(batches: list[dict[str, str]]) -> list[list[str]]:
    rows = [["runID", "machineID", "note", "index_hopping_rate", "method"]]
    for batch in batches:
        path = result_file(Path(batch["dir"]), "lambdaSD/index_hopping_rate.xls", "lambdaSD/results/index_hopping_rate.xls")
        if not path.exists():
            continue
        header, data = read_tsv(path, header=True)
        if "Index_Hopping_Rate_Percent" not in header or not data:
            raise ValidationError(f"Invalid lambdaSD rate table: {path}")
        rate = data[0][header.index("Index_Hopping_Rate_Percent")]
        rows.append([batch["run_id"], batch["run_id"].split("_")[1] if "_" in batch["run_id"] else "", "", f"{float(rate):.5f}%", "filtered_diagonal_neighbor_pm3"])
    return rows


def lambda_read_counts(batches: list[dict[str, str]]) -> list[list[str]]:
    rows = [["batch", "runID", "sampleID", "total_reads", "correctly_mapped", "misassigned_counts"]]
    for batch in batches:
        path = result_file(Path(batch["dir"]), "lambdaSD/read_align_count_filter.xls", "lambdaSD/results/read_align_count_filter.xls")
        if not path.exists():
            continue
        header, data = read_tsv(path, header=True)
        required = ["sampleID", "Total_Reads", "Correctly_Mapped", "Misassigned_Counts"]
        if any(name not in header for name in required):
            raise ValidationError(f"Invalid lambdaSD read-count table: {path}")
        idx = {name: header.index(name) for name in required}
        rows.extend([[batch["batch"], batch["run_id"], row[idx["sampleID"]], row[idx["Total_Reads"]], row[idx["Correctly_Mapped"]], row[idx["Misassigned_Counts"]]] for row in data])
    return rows


def mngs_qc(batches: list[dict[str, str]]) -> list[list[str]]:
    rows = [["batch", "runID", "library", "raw_reads", "raw_q20", "raw_q30", "raw_gc", "adapter_rate", "host_ratio", "dup_rate", "classified_reads", "classified_ratio", "unclassified_reads", "unclassified_ratio"]]
    for batch in batches:
        path = Path(batch["dir"]) / "mNGS/stat.txt"
        if not path.exists():
            continue
        _, data = read_tsv(path)
        for row in data:
            if len(row) < 14:
                raise ValidationError(f"Invalid mNGS QC row in {path}: expected 14 columns")
            rows.append([batch["batch"], batch["run_id"], *row[2:14]])
    return rows


def mngs_abundance(batches: list[dict[str, str]]) -> dict[str, list[list[str]]]:
    tables: dict[str, list[list[str]]] = {}
    for batch in batches:
        folder = Path(batch["dir"]) / "mNGS"
        if not folder.is_dir():
            continue
        for path in sorted(folder.glob("CN*.txt")):
            header, data = read_tsv(path, header=True)
            if len(header) < 2:
                raise ValidationError(f"Invalid mNGS abundance header: {path}")
            library = path.stem
            table = tables.setdefault(library, [["batch", "runID", "species", "taxid", "measured_pct", "expected_pct"]])
            for row in data:
                if len(row) < 2:
                    continue
                species = row[0]
                taxid = row[1] if len(row) >= 3 and row[1].strip().isdigit() else ""
                measured = row[-1] if len(row) == 3 else row[1]
                expected = row[-1] if len(row) > 3 else (row[2] if len(row) == 3 and not taxid else "")
                if len(row) == 3 and taxid:
                    measured, expected = row[2], ""
                table.append([batch["batch"], batch["run_id"], species, taxid, measured, expected])
    return tables


def tb2_qc(batches: list[dict[str, str]]) -> list[list[str]]:
    rows = [["batch", "runID", "prefix", "raw_reads", "clean_reads", "effective_rate", "adapter", "q20", "q30", "map", "on_target", "average_depth", "uniformity", "tb_drug", "ntm_drug", "ntm_id", "coinfection", "internal_control", "cov_1x", "cov_20x", "cov_50x", "cov_100x", "cov_200x", "cov_500x"]]
    for batch in batches:
        path = result_file(Path(batch["dir"]), "TB2/result/QC.stat.xls", "TB2/QC.stat.xls")
        if not path.exists():
            continue
        _, data = read_tsv(path, header=True)
        for row in data:
            if len(row) < 23:
                raise ValidationError(f"Invalid TB2 QC row in {path}")
            rows.append([batch["batch"], *row[1:23]])
    return rows


def tb2_depth(batches: list[dict[str, str]], kind: str) -> list[list[str]]:
    suffix = "_amplicon.depth.xls" if kind == "amplicon" else "_hot_depth.xls"
    base = ["batch", "runID", "library", "depth"]
    key_headers = ["amp"] if kind == "amplicon" else ["chr", "start", "end", "mutation"]
    rows = [base[:-1] + key_headers + ["depth"]]
    for batch in batches:
        candidates = sorted(result_dir(Path(batch["dir"]), "TB2/result", "TB2").glob(f"*{suffix}"))
        if not candidates:
            continue
        if len(candidates) != 1:
            raise ValidationError(f"Expected one {suffix} file for {batch['batch']}, got {len(candidates)}")
        header, data = read_tsv(candidates[0], header=True)
        key_count = len(key_headers)
        if len(header) <= key_count:
            raise ValidationError(f"Invalid TB2 depth table: {candidates[0]}")
        libraries = header[key_count:]
        for row in data:
            if len(row) < key_count:
                continue
            for offset, library in enumerate(libraries, start=key_count):
                rows.append([batch["batch"], batch["run_id"], library, *row[:key_count], row[offset] if offset < len(row) else "0"])
    return rows


def require_batch_has_data(batch: dict[str, str]) -> None:
    root = Path(batch["dir"])
    required = [root / "mNGS/stat.txt", result_file(root, "TB2/result/QC.stat.xls", "TB2/QC.stat.xls"), result_file(root, "lambdaSD/index_hopping_rate.xls", "lambdaSD/results/index_hopping_rate.xls")]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise ValidationError(f"Batch {batch['batch']} is incomplete; missing: {'; '.join(missing)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", help="Batch ID to validate and integrate, e.g. BGISEQ_0707")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()

    batches = parse_manifest(args.manifest)
    selected = next((item for item in batches if item["batch"].casefold() == args.batch.casefold()), None)
    if not selected:
        raise ValidationError(f"Batch {args.batch!r} is absent from {args.manifest}")
    require_batch_has_data(selected)

    # CNV/PGTA/tNGS are generated by their established project scripts before
    # this step.  Their aggregate TSVs become immutable local master CSVs.
    source = PROJECT_ROOT / "bin"
    tables: dict[Path, list[list[str]]] = {
        Path("cnvseq/CNVseq.csv"): existing_global_table(source / "cnv/stat.txt"),
        Path("pgt/PGTA.csv"): existing_global_table(source / "pgta/stat.txt"),
        Path("lambsd/汇总表.csv"): lambda_summary(batches),
        Path("lambsd/readscount表.csv"): lambda_read_counts(batches),
        Path("tNGS/阴阳性及hopping-新.csv"): existing_global_table(source / "tngs/test/result_stat.txt"),
        Path("tNGS/质控汇总表.csv"): existing_global_table(source / "tngs/test/tngs_stat.txt"),
        Path("tNGS/index hopping矩阵表-新.csv"): existing_global_table(source / "tngs/test/ap_result.txt"),
        Path("mNGS/mNGS质控汇总.csv"): mngs_qc(batches),
        Path("TB2/数据质控.csv"): tb2_qc(batches),
        Path("TB2/原始扩增子深度表.csv"): tb2_depth(batches, "amplicon"),
        Path("TB2/突变热点深度表.csv"): tb2_depth(batches, "hot"),
    }
    for library, rows in mngs_abundance(batches).items():
        tables[Path("mNGS") / f"{library}.csv"] = rows

    # The expected 15 deliverables are explicit: one QC table plus the four
    # control-library abundance tables.  A different library is accepted as an
    # additional table, but a missing control library fails early.
    expected = {
        Path("mNGS/CN002341-S01-D01-L02-UDI1.csv"),
        Path("mNGS/CN002386-S01-D01-L01-UDI2.csv"),
        Path("mNGS/CN002387-S01-D01-L01-UDI3.csv"),
        Path("mNGS/CN002388-S01-D01-L01-UDI4.csv"),
    }
    missing = expected - set(tables)
    if missing:
        raise ValidationError("Missing mNGS abundance tables: " + ", ".join(map(str, sorted(missing))))

    report = {"batch": selected["batch"], "run_id": selected["run_id"], "generated_at": datetime.now(timezone.utc).isoformat(), "tables": {}}
    if args.validate_only:
        for name, rows in tables.items():
            report["tables"][str(name)] = {"rows": len(rows), "columns": max(len(row) for row in rows)}
    else:
        for name, rows in tables.items():
            report["tables"][str(name)] = write_csv_atomic(args.out / name, rows)
        manifest = args.out / "runs" / f"{selected['batch']}.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        temp = manifest.with_suffix(".json.tmp")
        temp.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temp, manifest)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"VALIDATION FAILED: {exc}", file=sys.stderr)
        raise SystemExit(2)
