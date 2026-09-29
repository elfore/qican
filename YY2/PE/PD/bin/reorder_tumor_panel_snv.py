#!/usr/bin/env python3
"""
Reorder tumor-panel SNV rows by the fixed Feishu coordinate reference. 
cd /mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/SKII-JBJC-YY2/tumor_panel/result/SKII_3-CN001872-S02-D01-L02-UDI70
python3 /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/PE/PD/bin/reorder_tumor_panel_snv.py SKII_3-CN001872-S02-D01-L02-UDI70.snv_positive.txt SKII_3-CN001872-S02-D01-L02-UDI70.snv_positive_order.txt
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path


REFERENCE_FILES = (
    Path(__file__).with_name("SNV_all.csv"),
    Path(__file__).with_name("SNV突变频率结果.csv"),
)


@dataclass(frozen=True)
class Coordinate:
    order: int
    chromosome: str
    position: str


def coordinate_key(chromosome: str, position: str) -> tuple[str, str]:
    """Match chr1 and 1, while preserving the reference spelling on output."""
    chrom = chromosome.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.lower()

    pos = position.strip()
    if pos.isdigit():
        pos = str(int(pos))
    return chrom, pos


def open_text(path: Path):
    return path.open("r", encoding="utf-8-sig", newline="")


def detect_delimiter(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    with open_text(path) as handle:
        first_line = handle.readline()
    return "\t" if "\t" in first_line and "," not in first_line else ","


def read_reference(
    path: Path,
    delimiter: str,
    has_header: bool,
) -> dict[tuple[str, str], Coordinate]:
    reference: dict[tuple[str, str], Coordinate] = {}
    with open_text(path) as handle:
        rows = csv.reader(handle, delimiter=delimiter)
        if has_header:
            try:
                next(rows)
            except StopIteration as exc:
                raise ValueError(f"Reference file is empty: {path}") from exc

        for line_number, row in enumerate(rows, start=2 if has_header else 1):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) < 2:
                raise ValueError(
                    f"Reference row {line_number} has fewer than two columns: {row!r}"
                )
            chromosome, position = row[0].strip(), row[1].strip()
            if not chromosome or not position:
                raise ValueError(
                    f"Reference row {line_number} has an empty coordinate: {row!r}"
                )
            key = coordinate_key(chromosome, position)
            if key in reference:
                previous = reference[key]
                raise ValueError(
                    "Duplicate coordinate in reference: "
                    f"{chromosome}:{position} (orders "
                    f"{previous.order + 1} and {len(reference) + 1})"
                )
            reference[key] = Coordinate(
                order=len(reference),
                chromosome=chromosome,
                position=position,
            )
    if not reference:
        raise ValueError(f"Reference file has no coordinate rows: {path}")
    return reference


def read_result(
    path: Path,
    delimiter: str,
    has_header: bool,
) -> tuple[list[str] | None, list[list[str]]]:
    header: list[str] | None = None
    rows: list[list[str]] = []
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        for line_number, row in enumerate(reader, start=1):
            if not row or not any(cell.strip() for cell in row):
                continue
            if has_header and header is None:
                header = row
                continue
            if len(row) < 2:
                raise ValueError(
                    f"Result row {line_number} has fewer than two columns: {row!r}"
                )
            rows.append(row)
    if not rows:
        raise ValueError(f"Result file has no data rows: {path}")
    return header, rows


def reorder_rows(
    reference: dict[tuple[str, str], Coordinate],
    rows: list[list[str]],
    input_name: str,
) -> tuple[list[list[str]], int]:
    indexed_rows: list[tuple[int, int, list[str]]] = []
    seen: set[tuple[str, str]] = set()
    unknown: list[str] = []
    duplicate: list[str] = []

    for row_number, row in enumerate(rows, start=1):
        key = coordinate_key(row[0], row[1])
        coordinate = reference.get(key)
        if coordinate is None:
            unknown.append(f"{row[0]}:{row[1]} (result row {row_number})")
            continue
        if key in seen:
            duplicate.append(f"{row[0]}:{row[1]} (result row {row_number})")
            continue
        seen.add(key)

        normalized = list(row)
        normalized[0] = coordinate.chromosome
        normalized[1] = coordinate.position
        indexed_rows.append((coordinate.order, row_number, normalized))

    if unknown:
        preview = ", ".join(unknown[:10])
        suffix = " ..." if len(unknown) > 10 else ""
        raise ValueError(
            f"{input_name}: {len(unknown)} result coordinate(s) are absent from "
            f"the reference: {preview}{suffix}"
        )
    if duplicate:
        preview = ", ".join(duplicate[:10])
        suffix = " ..." if len(duplicate) > 10 else ""
        raise ValueError(
            f"{input_name}: duplicate result coordinate(s): {preview}{suffix}"
        )

    missing_count = len(reference) - len(seen)
    indexed_rows.sort(key=lambda item: (item[0], item[1]))
    return [row for _, _, row in indexed_rows], missing_count


def expand_to_reference(
    reference: dict[tuple[str, str], Coordinate],
    rows: list[list[str]],
) -> tuple[list[list[str]], int]:
    """Emit every reference coordinate and blank fields for absent result rows."""
    result_by_key: dict[tuple[str, str], list[str]] = {}
    for row in rows:
        key = coordinate_key(row[0], row[1])
        coordinate = reference[key]
        normalized = list(row)
        normalized[0] = coordinate.chromosome
        normalized[1] = coordinate.position
        result_by_key[key] = normalized

    width = max((len(row) for row in rows), default=2)
    output: list[list[str]] = []
    missing_count = 0
    for key, coordinate in sorted(
        reference.items(), key=lambda item: item[1].order
    ):
        row = result_by_key.get(key)
        if row is None:
            missing_count += 1
            row = [coordinate.chromosome, coordinate.position] + [""] * (width - 2)
        output.append(row)
    return output, missing_count


def write_result(
    output: Path,
    header: list[str] | None,
    rows: list[list[str]],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_MINIMAL,
        )
        if header is not None:
            writer.writerow(header)
        writer.writerows(rows)


def resolve_reference_file() -> Path:
    for path in REFERENCE_FILES:
        if path.is_file():
            return path
    names = ", ".join(path.name for path in REFERENCE_FILES)
    raise FileNotFoundError(
        f"Reference file not found beside script; expected one of: {names}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reorder a tumor-panel SNV result by the fixed Feishu coordinate "
            "reference kept beside this script."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Input SNV result file.",
    )
    parser.add_argument(
        "output",
        type=Path,
        help="Output file.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"Input file does not exist: {args.input}")

    try:
        reference_file = resolve_reference_file()
        reference_delimiter = detect_delimiter(reference_file, "auto")
        reference = read_reference(
            reference_file,
            reference_delimiter,
            has_header=True,
        )
        header, rows = read_result(
            args.input,
            "\t",
            has_header=False,
        )
        validated_rows, _ = reorder_rows(
            reference,
            rows,
            str(args.input),
        )
        reordered, missing_count = expand_to_reference(reference, validated_rows)
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    try:
        write_result(
            args.output,
            header,
            reordered,
        )
    except OSError as exc:
        print(f"ERROR: failed to write {args.output}: {exc}", file=sys.stderr)
        return 1

    print(f"reference_rows={len(reference)}")
    print(f"result_rows={len(rows)}")
    print(f"missing_reference_rows={missing_count}")
    print(f"output_rows={len(reordered)}")
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
