#!/usr/bin/env python3
"""YY1 controlled one-click aggregation workflow.

The workflow is deliberately two-phase.  It always builds and validates the
15 local master CSV files first.  Feishu writes are opt-in via
``--apply-feishu`` so that an incomplete batch can never clear an online table.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


YY1 = Path("/mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY1")


def run(argv: list[str], *, cwd: Path | None = None) -> None:
    print("[RUN]", " ".join(argv), flush=True)
    subprocess.run(argv, cwd=cwd, check=True)


def refresh_small_aggregate_tables() -> None:
    """Rebuild lightweight aggregate TSVs from completed result directories."""
    run(["python", "qc_stat.py", "-infile", "../YY1_info.txt", "-infile2", "cnv.txt", "-outfile", "stat.txt"], cwd=YY1 / "bin/cnv")
    run(["python", "qc_stat.py", "-infile", "../YY1_info.txt", "-infile2", "cnv.txt", "-outfile", "stat.txt"], cwd=YY1 / "bin/pgta")
    tngs = YY1 / "bin/tngs/test"
    run(["python", "ap_stat.py", "-infile", "../../YY1_info.txt", "-pos", "info.txt", "-outfile", "ap_result.txt", "-pos2", "sp_ap_info.txt"], cwd=tngs)
    run(["python", "result_check.py", "-infile", "../../YY1_info.txt", "-pos", "result_info_new.txt", "-outfile", "result_stat.txt"], cwd=tngs)
    run(["python", "qc_stat.py", "-infile", "../../YY1_info.txt", "-pos", "result_info_new.txt"], cwd=tngs)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", help="Batch ID registered in bin/YY1_info.txt")
    parser.add_argument("--refresh-aggregate-tsv", action="store_true", help="Run the existing lightweight CNV/PGTA/tNGS aggregate scripts before building CSVs.")
    parser.add_argument("--apply-feishu", action="store_true", help="After a successful CSV build, update all 15 Feishu sheets and refresh XLSX backups.")
    args = parser.parse_args()

    if args.refresh_aggregate_tsv:
        refresh_small_aggregate_tables()

    # No online write precedes this command.  It creates atomically replaced
    # local CSV files plus a run manifest, or fails without starting sync.
    run(["python", str(YY1 / "bin/build_assessment_csv.py"), args.batch])

    if not args.apply_feishu:
        print("[OK] Local CSV master tables are validated. Feishu unchanged (pass --apply-feishu to sync).")
        return 0

    # Compatibility bridge for the original 11 sheets.  It remains isolated so
    # the legacy full-rebuild behaviour can be retired sheet-by-sheet once each
    # master CSV projection has been migrated.
    run(["python", str(YY1 / "bin/rebuild_and_sync.py")])
    run(["python", str(YY1 / "bin/sync_mngs_abundance.py"), args.batch, "--apply"])
    run(["python", str(YY1 / "bin/export_backup_with_retry.py")])
    print("[OK] Feishu sync and six XLSX backups completed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"[ERROR] Command failed with exit code {exc.returncode}: {exc.cmd}", file=sys.stderr)
        raise SystemExit(exc.returncode)
