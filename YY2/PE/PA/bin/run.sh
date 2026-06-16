#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="/mnt/gpfs1/Users/yangjinxurong/projects/qican/YY2/PE/PA/bin"
KB_FILE="/mnt/gpfs1/Users/caiyilun/test/qican_stat/chd/kownledge_base.txt"

cd "$SCRIPT_DIR"

python stat.py -infile ../../YY2_info.txt -pos "${KB_FILE}" -outfile ../qc_stat.txt
python sample_stat.py -infile ../../YY2_info.txt -pos "${KB_FILE}" -outfile ../vars_stat.txt
python sample_stat_allv.py -infile ../../YY2_info.txt -pos "${KB_FILE}" -outfile ../vars_stat_all.txt
python stat_dp_1.py -infile ../../YY2_info_dp.txt -pos "${KB_FILE}" -outfile ../dp_stat.txt
