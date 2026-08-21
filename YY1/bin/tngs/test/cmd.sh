python ap_stat.py -infile ../../YY1_info.txt -pos info.txt -outfile ap_result.txt -pos2 sp_ap_info.txt
python result_check.py -infile ../../YY1_info.txt -pos result_info_new.txt -outfile result_stat.txt
python qc_stat.py -infile ../../YY1_info.txt -pos result_info_new.txt 