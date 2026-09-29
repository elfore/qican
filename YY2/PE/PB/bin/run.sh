#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="/mnt/gpfs1/Users/yangjinxurong/pipeline/qican/YY2/PE/PB/bin/"
cd "$SCRIPT_DIR"

python stat.py -infile ../../YY2_info.txt -outfile ../result.txt
python sample_stat.py -infile ../../YY2_info.txt -pos /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/YY2-PL24-022-01/PB/result/combine_mut/PL2402201_LC002901-S01-D01-L12-C208C160_combine_rs.vcf -outfile ../vars_stat.txt
python sample_stat_cyp_ty.py -infile ../../YY2_info.txt -outfile ../stat_cyp_ty.txt
python sample_stat_hla_new.py -infile ../../YY2_info.txt -outfile ../stat_hla.txt
python sample_stat_allv.py -infile ../../YY2_info.txt -pos /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/YY2-PL24-022-01/PB/result/combine_mut/PL2402201_LC002901-S01-D01-L12-C208C160_combine_rs.vcf -outfile ../vars_stat_all.txt -infile2 HLA.txt
python sample_stat_hla_ty.py -infile ../../YY2_info.txt -outfile ../stat_hla_ty.txt -infile2 HLA.txt
python stat_dp_1.py -infile ../../YY2_info_dp.txt -pos /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/YY2-PL24-022-01/PB/result/combine_mut/PL2402201_LC002901-S01-D01-L12-C208C160_combine_rs.vcf -outfile ../dp_stat.txt
