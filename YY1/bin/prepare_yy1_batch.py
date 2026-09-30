#!/usr/bin/env python3
"""
Create a YY1 batch workspace from the cpScript template.
/mnt/gpfs1/Users/yangjinxurong/software/miniconda3/envs/singlecell/bin/python3 prepare_yy1_batch.py --batchid 260920142132_B174_SKII-JBJC-YY1-260920142133
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

OUTPUT_ROOT = Path("/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/SE75")

YY1_ROOT = Path("/mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1")
TEMPLATE_ROOT = YY1_ROOT / "cpScript"

# These are the placeholders used by the existing YY1 template scripts.
BATCH_PLACEHOLDERS = ("batchID", "BATCHID")
SKID_PLACEHOLDERS = ("skID", "skid", "SKID")
FASTQ_RE = re.compile(r"""(?P<path>[^\s"'`<>|;&()]+?\.(?:fq|fastq)\.gz)""")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy YY1 cpScript templates, replace batch/skid placeholders, and check FASTQ paths."
    )
    parser.add_argument("--batchid", required=True, help="Full batch ID.")
    parser.add_argument(
        "--skid",
        help="SKID. If omitted, parse it from <prefix>_<instrument>_<skid>-<suffix>.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
        help=f"YY1 result root (default: {OUTPUT_ROOT}).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow writing into an existing result directory.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report planned actions without copying or changing files.",
    )
    return parser.parse_args()


def infer_skid(batchid: str) -> str:
    # Greedy middle group preserves SKIDs that themselves contain underscores.
    match = re.match(r"^[^_]+_[^_]+_(.+)-[^-]+$", batchid)
    if not match:
        raise ValueError(
            "Cannot infer skid from batchid. Please provide it explicitly with --skid."
        )
    skid = match.group(1)
    if not skid:
        raise ValueError("The inferred skid is empty.")
    return skid


def replace_placeholders(text: str, batchid: str, skid: str) -> str:
    for placeholder in BATCH_PLACEHOLDERS:
        text = text.replace(placeholder, batchid)
    for placeholder in SKID_PLACEHOLDERS:
        text = text.replace(placeholder, skid)
    return text


def extract_fastq_paths(text: str) -> list[str]:
    paths: list[str] = []
    for match in FASTQ_RE.finditer(text):
        path = match.group("path").rstrip(".,:")
        if path not in paths:
            paths.append(path)
    return paths


def load_sample_file(sample_file: Path) -> list[str]:
    """Read non-empty, non-comment lines from a sample list file."""
    if not sample_file.is_file():
        return []
    try:
        lines = sample_file.read_text(encoding="utf-8").splitlines()
        return [
            line.strip()
            for line in lines
            if line.strip() and not line.strip().startswith("#")
        ]
    except (UnicodeDecodeError, OSError):
        return []


def resolve_fastq_candidates(file_path: Path, raw_path: str) -> list[str]:
    """
    Resolve FASTQ paths containing variables like ${i}, $i, ${sample}, $sample
    by substituting sample names from the 'sample' file in the same directory.
    """
    if any(var in raw_path for var in ("${i}", "$i", "${sample}", "$sample")):
        sample_file = file_path.parent / "sample"
        if not sample_file.is_file():
            for alt in ("sample.list", "sample.txt"):
                alt_path = file_path.parent / alt
                if alt_path.is_file():
                    sample_file = alt_path
                    break
        samples = load_sample_file(sample_file)
        if samples:
            resolved = []
            for s in samples:
                item = (
                    raw_path.replace("${i}", s)
                    .replace("$i", s)
                    .replace("${sample}", s)
                    .replace("$sample", s)
                )
                resolved.append(item)
            return resolved
    return [raw_path]


def check_fastq_paths(result_root: Path) -> tuple[list[str], list[str]]:
    missing: list[str] = []
    unresolved: list[str] = []
    for path in sorted(result_root.rglob("*")):
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for raw_fastq in extract_fastq_paths(text):
            candidates = resolve_fastq_candidates(path, raw_fastq)
            for fastq_path in candidates:
                if "$" in fastq_path or "{" in fastq_path or "}" in fastq_path:
                    unresolved.append(f"{path}: {fastq_path}")
                    continue
                candidate = Path(fastq_path)
                if not candidate.is_absolute():
                    candidate = path.parent / candidate
                if not candidate.exists():
                    missing.append(str(candidate))
    return sorted(set(missing)), sorted(set(unresolved))


def copy_template(destination: Path) -> None:
    """Copy the cpScript tree contents into destination."""
    destination.mkdir(parents=True, exist_ok=True)
    # copytree on the root copies its contents into an existing destination,
    # including hidden entries
    # 临时禁用 copystat，阻止 copytree 覆写目录的历史时间戳
    orig_copystat = shutil.copystat
    try:
        shutil.copystat = lambda src, dst, **kw: None
        shutil.copytree(
            TEMPLATE_ROOT,
            destination,
            dirs_exist_ok=True,
            copy_function=shutil.copy,  # 用 copy 替代 copy2，不拷贝文件的历史时间
        )
    finally:
        shutil.copystat = orig_copystat


def replace_in_files(result_root: Path, batchid: str, skid: str) -> int:
    changed = 0
    for path in result_root.rglob("*"):
        if not path.is_file():
            continue
        try:
            original = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        updated = replace_placeholders(original, batchid, skid)
        if updated != original:
            with path.open("w", encoding="utf-8", newline="") as f:
                f.write(updated)
            changed += 1
    return changed


def main() -> int:
    args = parse_args()
    batchid = args.batchid
    skid = args.skid or infer_skid(batchid)
    output_root = args.output_root.resolve()
    result_root = output_root / skid

    if not re.fullmatch(r"[A-Za-z0-9._-]+", skid):
        raise ValueError(f"Unsafe skid for directory name: {skid!r}")
    if not TEMPLATE_ROOT.is_dir():
        raise FileNotFoundError(f"Template directory does not exist: {TEMPLATE_ROOT}")
    if result_root.exists() and not args.force:
        # An existing directory is allowed when it is empty; this avoids accidental
        # overwrites while still permitting a normal first run after mkdir -p.
        if any(result_root.iterdir()):
            raise FileExistsError(
                f"Result directory is not empty: {result_root}. Use --force to merge-copy."
            )

    print(f"[INFO] batchid={batchid}")
    print(f"[INFO] skid={skid}")
    print(f"[INFO] result={result_root}")
    print(f"[INFO] template={TEMPLATE_ROOT}")
    if args.dry_run:
        print("[DRY-RUN] No files will be copied or modified.")
        return 0

    result_root.mkdir(parents=True, exist_ok=True)
    copy_template(result_root)
    changed = replace_in_files(result_root, batchid, skid)
    missing, unresolved = check_fastq_paths(result_root)

    print(f"[INFO] replaced files={changed}")
    if missing:
        print("[MISSING_FASTQ]")
        for path in missing:
            print(path)
    else:
        print("[INFO] No missing literal FASTQ paths found.")
    if unresolved:
        print("[UNRESOLVED_FASTQ]")
        for item in unresolved:
            print(item)
    print(f"[SUMMARY] missing_fastq={len(missing)} unresolved_fastq={len(unresolved)}")
    return 1 if missing else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(2)
