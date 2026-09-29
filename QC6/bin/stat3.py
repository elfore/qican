import argparse
from calendar import c
import os
import sys
import re
import pandas as pd
from pandas import DataFrame

# cd ../results
# /mnt/gpfs/Users/yangjinxurong/software/miniconda3/envs/stats/bin/python /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/QC6/bin/stat3.py -indir /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/05.Qican6/dingzhi_final -infile /mnt/gpfs1/Users/yangjinxurong/pipeline/qican/QC6/bin/aa.txt

def to_number(x: str):
    s = x.strip()
    if s.endswith('%'):
        s_no_percent = s[:-1].strip()  # 去掉 '%'
        try:
            # 要除以 100，这样 "98.54%" -> 0.9854
            return float(s_no_percent) / 100.0
        except ValueError:
            return x

    # 如果不是百分号，再判断是否是整数/浮点/科学计数法
    # 匹配整数: 123 或 -123 等
    if re.fullmatch(r"[+-]?\d+", s):
        return int(s)
    # 匹配浮点: 123.456, -1.2, 3.5e10, -2E-2 等
    elif re.fullmatch(r"[+-]?\d+(\.\d+)?([eE][+-]?\d+)?", s):
        return float(s)
    else:
        return x


def parse_line_as_numbers(line: str):
    """
    将一行（以 \t 分割）中的各字段尽量转换为数值型。
    """
    fields = line.strip("\n").split("\t")
    fields = [to_number(field) for field in fields]
    return fields


def add_min_max_rows(df: pd.DataFrame) -> pd.DataFrame:
    numeric_cols = df.select_dtypes(include=[int, float]).columns

    min_row = []
    max_row = []
    for index,col in enumerate(df.columns):
        if index==0:
            min_row.append("Min")
            max_row.append("Max")
        else:
            if col in numeric_cols:
                col_min = df[col].min()
                col_max = df[col].max()
                min_row.append(col_min)
                max_row.append(col_max)
            else:
                min_row.append('-')
                max_row.append('-')
    
    # 在 df 最后添加两行（并给它们设定自定义索引）
    df.loc["Min"] = min_row
    df.loc["Max"] = max_row
    return df

def contains_all_elements(text, elements_list):
    return all(element in text for element in elements_list)


def stat(infile, indir, outdir):
    head1 = [
        "run_id","version","platform","species","raw_reads","raw_bases","raw_len_read1",
        "raw_len_read2","raw_q20","raw_q30","raw_q40","raw_gc_content","fq_duprate",
        "reads_with_adapter_rate","clean_reads","clean_bases","clean_len_read1","clean_len_read2",
        "clean_q20","clean_q30","clean_q40","clean_gc_content","Mapping_rate","Mapping_rate_MQ30",
        "properly_paired_rate","Duplication_Rate_bam","average_depth_mq30","uniformity_mq30_20%",
        "uniformity_mq30_50%","Cov_at_1X_mq30","Cov_at_20X_mq30","Cov_at_50X_mq30",
        "Cov_at_100X_mq30","lowGC-normdepth-bases","evenGC-normdepth-bases","highGC-normdepth-bases",
        "depth-cv-bases","lowGC-normdepth-reads","evenGC-normdepth-reads","highGC-normdepth-reads",
        "depth-cv-reads","lowGC-diff-reads","evenGC-diff-reads","highGC-diff-reads",
        "clean_mean_tiletup_1000","clean_median_tiletup_1000","clean_tiledup_cv_1000","Mismatch_rate",
        "1bp_ins_rate","multibase_ins_rate","total_ins_rate","1bp_del_rate","multibase_del_rate",
        "total_del_rate","Total_error_rates","totalpos",">0.1%_error_pos_rate_all",">1%_error_pos_rate_all",
        ">5%_error_pos_rate_all",">10%_error_pos_rate_all",">0.1%_error_pos_rate_SNV",
        ">1%_error_pos_rate_SNV",">5%_error_pos_rate_SNV",">10%_error_pos_rate_SNV","totalpos_rm_dpl_str",
        ">0.1%_error_pos_rate_all_rm_dpl_str",">1%_error_pos_rate_all_rm_dpl_str",">5%_error_pos_rate_all_rm_dpl_str",
        ">10%_error_pos_rate_all_rm_dpl_str",">0.1%_error_pos_rate_SNV_rm_dpl_str",">1%_error_pos_rate_SNV_rm_dpl_str",
        ">5%_error_pos_rate_SNV_rm_dpl_str",">10%_error_pos_rate_SNV_rm_dpl_str","totalpos_dpl_str",
        ">0.1%_error_pos_rate_all_dpl_str",">1%_error_pos_rate_all_dpl_str",">5%_error_pos_rate_all_dpl_str",
        ">10%_error_pos_rate_all_dpl_str",">0.1%_error_pos_rate_SNV_dpl_str",">1%_error_pos_rate_SNV_dpl_str",
        ">5%_error_pos_rate_SNV_dpl_str",">10%_error_pos_rate_SNV_dpl_str","sensitivity_SNV","precision_SNV",
        "F1score_SNV","sensitivity_indel","precision_indel","F1score_indel","total_fp","snv_fp(10%)",
        "1bp_ins_fp(10%)","1bp_del_fp(10%)","1bp_indel_fp(10%)","snv_fp(5%)","1bp_ins_fp(5%)",
        "1bp_del_fp(5%)","1bp_indel_fp(5%)","snv_fp(1%)","1bp_ins_fp(1%)","1bp_del_fp(1%)","1bp_indel_fp(1%)",
        "snv_fp(0.5%)","1bp_ins_fp(0.5%)","1bp_del_fp(0.5%)","1bp_indel_fp(0.5%)","Pathogenic",
        "Pathogenic_1bpindel","most_common","mean_insertsize","median_insertsize","insertsize_std","Skewness"
    ]
    head2 = [
        "run_id","version","platform","species","Ratio_relative_to_350bp","raw_reads","raw_bases",
        "raw_len_read1","raw_len_read2","raw_q20","raw_q30","raw_q40","raw_gc_content","fq_duprate",
        "reads_with_adapter_rate","clean_reads","clean_bases","clean_len_read1","clean_len_read2",
        "clean_q20","clean_q30","clean_q40","clean_gc_content","Mapping_rate","Mapping_rate_MQ30",
        "properly_paired_rate","Duplication_Rate_bam","average_depth_mq30","uniformity_mq30_20%",
        "uniformity_mq30_50%","Cov_at_1X_mq30","Cov_at_20X_mq30","Cov_at_50X_mq30","Cov_at_100X_mq30",
        "lowGC-normdepth-bases","evenGC-normdepth-bases","highGC-normdepth-bases","depth-cv-bases",
        "lowGC-normdepth-reads","evenGC-normdepth-reads","highGC-normdepth-reads","depth-cv-reads",
        "lowGC-diff-reads","evenGC-diff-reads","highGC-diff-reads","clean_mean_tiletup_1000",
        "clean_median_tiletup_1000","clean_tiledup_cv_1000","Mismatch_rate","1bp_ins_rate",
        "multibase_ins_rate","total_ins_rate","1bp_del_rate","multibase_del_rate","total_del_rate",
        "Total_error_rates","totalpos",">0.1%_error_pos_rate_all",">1%_error_pos_rate_all",">5%_error_pos_rate_all",
        ">10%_error_pos_rate_all",">0.1%_error_pos_rate_SNV",">1%_error_pos_rate_SNV",">5%_error_pos_rate_SNV",
        ">10%_error_pos_rate_SNV","totalpos_rm_dpl_str",">0.1%_error_pos_rate_all_rm_dpl_str",
        ">1%_error_pos_rate_all_rm_dpl_str",">5%_error_pos_rate_all_rm_dpl_str",">10%_error_pos_rate_all_rm_dpl_str",
        ">0.1%_error_pos_rate_SNV_rm_dpl_str",">1%_error_pos_rate_SNV_rm_dpl_str",
        ">5%_error_pos_rate_SNV_rm_dpl_str",">10%_error_pos_rate_SNV_rm_dpl_str","totalpos_dpl_str",
        ">0.1%_error_pos_rate_all_dpl_str",">1%_error_pos_rate_all_dpl_str",">5%_error_pos_rate_all_dpl_str",
        ">10%_error_pos_rate_all_dpl_str",">0.1%_error_pos_rate_SNV_dpl_str",">1%_error_pos_rate_SNV_dpl_str",
        ">5%_error_pos_rate_SNV_dpl_str",">10%_error_pos_rate_SNV_dpl_str","sensitivity_SNV","precision_SNV",
        "F1score_SNV","sensitivity_indel","precision_indel","F1score_indel","total_fp","snv_fp(10%)",
        "1bp_ins_fp(10%)","1bp_del_fp(10%)","1bp_indel_fp(10%)","snv_fp(5%)","1bp_ins_fp(5%)",
        "1bp_del_fp(5%)","1bp_indel_fp(5%)","snv_fp(1%)","1bp_ins_fp(1%)","1bp_del_fp(1%)","1bp_indel_fp(1%)",
        "snv_fp(0.5%)","1bp_ins_fp(0.5%)","1bp_del_fp(0.5%)","1bp_indel_fp(0.5%)","Pathogenic",
        "Pathogenic_1bpindel","most_common","mean_insertsize","median_insertsize","insertsize_std","Skewness"
    ]
    head3 = [
        "runID","Version","Platform","type","L0_rawreads_ratio","L100_rawreads_ratio","L200_rawreads_ratio",
        "L300_rawreads_ratio","L400_rawreads_ratio","L500_rawreads_ratio","L600_rawreads_ratio",
        "L700_rawreads_ratio","L800_rawreads_ratio","L1000_rawreads_ratio","L1200_rawreads_ratio",
        "L0_mapped_ratio","L100_mapped_ratio","L200_mapped_ratio","L300_mapped_ratio","L400_mapped_ratio",
        "L500_mapped_ratio","L600_mapped_ratio","L700_mapped_ratio","L800_mapped_ratio","L1000_mapped_ratio",
        "L1200_mapped_ratio","L0_inertsize0_ratio","L100_inertsize0_ratio","L200_inertsize0_ratio",
        "L300_inertsize0_ratio","L400_inertsize0_ratio","L500_inertsize0_ratio","L600_inertsize0_ratio",
        "L700_inertsize0_ratio","L800_inertsize0_ratio","L1000_inertsize0_ratio","L1200_inertsize0_ratio"
    ]
    head4 = [
        "runid","version","platform","type","PA227","PA364","PA384","weiyan_TB1","weiyan_TB2",
        "PA227_fp_hotspot_count","PA364_fp_hotspot_count","PA384_fp_hotspot_count","PA364_11210898G>C",
        "PA227_mut_count_gt001","PA227_mut_count_gt005","PA227_mut_count_gt010",
        "PA364_mut_count_gt001","PA364_mut_count_gt005","PA364_mut_count_gt010","PA384_mut_count_gt001",
        "PA384_mut_count_gt005","PA384_mut_count_gt010","weiyan_TB1_mut_count_gt001","weiyan_TB1_mut_count_gt005",
        "weiyan_TB1_mut_count_gt010","weiyan_TB2_mut_count_gt001","weiyan_TB2_mut_count_gt005","weiyan_TB2_mut_count_gt010"
    ]
    head5 = [
        "runid","version","platform","type",">=80%_cov_amplicon_ratio",">=90%_cov_amplicon_ratio",
        ">=95%_cov_amplicon_ratio","KM_ac_10_0","KM_ac_10_10","KM_ac_10_2","KM_ac_10_5",
        "KM_ac_12_0","KM_ac_12_12","KM_ac_14_0","KM_ac_14_14","KM_ac_16_0","KM_ac_16_16","KM_ac_16_8",
        "KM_ac_18_0","KM_ac_20_2","KM_ca_10_0","KM_ca_10_2","KM_gt_10_0","KM_gt_10_2","KM_gt_10_5",
        "KM_gt_14_14","KM_gt_16_2","KM_tg_10_0_2","KM_tg_10_0_3","KM_tg_10_2","KM_tg_12_0","KM_tg_14_0",
        "KM_tg_16_0","KM_tg_16_2","KM_tg_18_0","KM_tg_20_0","KM_tg_20_2","SW_at_10_10","SW_at_14_0",
        "SW_at_16_0","SW_at_18_0","SW_at_20_0","SW_cg_10_0","SW_cg_10_2","SW_cg_10_5","SW_gc_10_2",
        "SW_gc_10_5","SW_ta_10_0","SW_ta_12_0","SW_ta_16_0","SW_ta_16_2","SW_ta_18_0","SW_ta_20_0",
        "YR_ag_12_0","YR_ag_14_0","YR_ag_16_0","YR_ag_16_2","YR_ag_16_8","YR_ag_18_0","YR_ag_20_2",
        "YR_ct_10_10","YR_ct_10_5","YR_ct_14_14","YR_ct_16_8","YR_ga_10_0","YR_ga_10_2","YR_ga_12_0",
        "YR_tc_10_0","YR_tc_16_0","YR_tc_18_0"
    ]
    head6 = [
        "run_id","version","platform","species","raw_reads","raw_bases","raw_len_read1","raw_len_read2",
        "raw_q20","raw_q30","raw_q40","raw_gc_content","fq_duprate","reads_with_adapter_rate",
        "clean_reads","clean_bases","clean_len_read1","clean_len_read2","clean_q20","clean_q30","clean_q40",
        "clean_gc_content","Mapping_rate","Mapping_rate_MQ30","properly_paired_rate","Duplication_Rate_bam",
        "average_depth_mq30","uniformity_mq30_20%","uniformity_mq30_50%","Cov_at_1X_mq30","Cov_at_20X_mq30",
        "Cov_at_50X_mq30","Cov_at_100X_mq30","Capture_Rate_on_Bases","Capture_Rate_on_Reads","FOLD_80_BASE_PENALTY",
        "AT_DROPOUT","GC_DROPOUT","Mismatch_rate","1bp_ins_rate","multibase_ins_rate","total_ins_rate",
        "1bp_del_rate","multibase_del_rate","total_del_rate","Total_error_rates","totalpos",">1%_error_pos_rate_all",
        ">1%_error_pos_rate_SNV","totalpos_rm_dpl_str",">1%_error_pos_rate_all_rm_dpl_str",
        ">1%_error_pos_rate_SNV_rm_dpl_str","totalpos_dpl_str",">1%_error_pos_rate_all_dpl_str",
        ">1%_error_pos_rate_SNV_dpl_str","lowGC-normdepth-bases","evenGC-normdepth-bases","highGC-normdepth-bases",
        "depth-cv-bases","lowGC-normdepth-reads","evenGC-normdepth-reads","highGC-normdepth-reads","depth-cv-reads",
        "lowGC-diff-reads","evenGC-diff-reads","highGC-diff-reads","clean_mean_tiletup_1000","clean_median_tiletup_1000",
        "clean_tiledup_cv_1000","sensitivity","specificity","precision","accuracy","f1_score","rmse","mae",
        "vaf1-5_pearson_corr","vaf5_10_pearson_corr","vaf10_100%","vaf1_5_spearman_corr","vaf5_10_spearman_corr",
        "vaf1-100%_spearman_corr","FP_SNP","FN_SNP","Precision_SNP","Sensitivity_SNP","F_measure_SNP","FP_INDEL",
        "FN_INDEL","Precision_INDEL","Sensitivity_INDEL","F_measure_INDEL","FP_SNP_0.5","FP_SNP_1","FP_SNP_10",
        "FP_INDEL_0.5","FP_INDEL_1","FP_INDEL_10","1bp_INDEL_FP_0.5","1bp_INDEL_FP_1","1bp_INDEL_FP_10","C-T","G-A",
        "G-T","C-A","T-C","A-G","A-C","A-T","C-G","G-C","T-A","T-G","most_common","mean_insertsize","median_insertsize",
        "insertsize_std","Skewness","mean-nordepth(dpl<=10)","mean-nordepth(10<dpl<=15)","mean-nordepth(15<dpl<=20)",
        "mean-nordepth(20<dpl<=25)","mean-nordepth(25<dpl<=30)","mean-nordepth(dpl>30)","mean-trimratio(dpl<=10)",
        "mean-trimratio(10<dpl<=15)","mean-trimratio(15<dpl<=20)","mean-trimratio(20<dpl<=25)",
        "mean-trimratio(25<dpl<=30)","mean-trimratio(dpl>30)","MSI_Max_depth","MSI_True_ratio",
        "MSI_True_ratio(in available)"
    ]

    with open(infile, "r") as IN:
        for line in IN:
            print(line)
            # 计数器
            num_350 = 0; num_tn5 = 0; num_small = 0; num_diff = 0
            num_burk=0; num_aur=0; num_brca=0; num_lambda=0
            num_len=0; num_spec=0; num_dpl=0; num_wes=0

            # 建立空的 DataFrame，后面按行追加
            df_350=DataFrame(columns = head1)
            df_tn5=DataFrame(columns = head1)
            df_small=DataFrame(columns = head1)
            df_diff=DataFrame(columns = head2)
            df_burk=DataFrame(columns = head1)
            df_aur=DataFrame(columns = head1)
            df_brca=DataFrame(columns = head1)
            df_lambda=DataFrame(columns = head1)
            df_len=DataFrame(columns = head3)
            df_spec=DataFrame(columns = head4)
            df_dpl=DataFrame(columns = head5)
            df_wes=DataFrame(columns = head6)

            seq_type=line.strip()
            writer = pd.ExcelWriter(f"{outdir}/{seq_type}_stat.xlsx")

            # 1) 读取 samples_qc_sum.txt
            list_qc = os.popen(f"ls {indir}/{seq_type}/*/*/result/QC_all/*samples_qc_sum.txt").read().strip("\n").split("\n")
            for i in list_qc:
                if not i:
                    continue
                if i.find("PL20250922-3")!=-1:continue  # 过滤掉这个样本
                with open(i, "r") as QC:
                    next(QC)  # 跳过表头
                    reads_350 = 0
                    for qc_line in QC:
                        lines = parse_line_as_numbers(qc_line)
                        if not lines:
                            continue
                        if len(lines) < 5:
                            continue
                        sample_name = lines[3]
                        # 根据 sample_name 判断分类
                        if contains_all_elements(sample_name,['Escherichia_coli','CN001913','UDI63']):
                            num_350 += 1
                            df_350.loc[num_350] = lines
                            reads_350 = lines[4]  # raw_reads
                        if contains_all_elements(sample_name,['Escherichia_coli','CN001913','N501N719']):
                            num_tn5 += 1
                            df_tn5.loc[num_tn5] = lines
                        if contains_all_elements(sample_name,['Escherichia_coli','CN001913','UDI13']):
                            num_small += 1
                            df_small.loc[num_small] = lines
                        suffixes = ['UDI63', 'UDI60', 'UDI61', 'UDI66']
                        if any(contains_all_elements(sample_name, ['Escherichia_coli', 'CN001913', suffix]) for suffix in suffixes):
                            num_diff += 1
                            if reads_350:
                                diff_nums = lines[4] / reads_350 if reads_350 != 0 else 0
                            else:
                                diff_nums = 0
                            # 在第 5 列插入 diff_nums
                            text_list = lines[0:4] + [diff_nums] + lines[4:]
                            df_diff.loc[num_diff] = text_list
                        if contains_all_elements(sample_name,['Burkholderia_multivorans','CN001914']):
                            num_burk += 1
                            df_burk.loc[num_burk] = lines
                        if contains_all_elements(sample_name,['Staphylococcus_aureus','CN001915']):
                            num_aur += 1
                            df_aur.loc[num_aur] = lines
                        if contains_all_elements(sample_name,['BRCA']):
                            num_brca += 1
                            df_brca.loc[num_brca] = lines
                        if contains_all_elements(sample_name,['lambda']):
                            num_lambda += 1

                            df_lambda.loc[num_lambda] = lines

            # 2) 读取 length_bias_stat_qc.xls
            list_len = os.popen(f"ls {indir}/{seq_type}/*/*/result/QC_all/*length_bias_stat_qc.xls").read().strip("\n").split("\n")
            for i in list_len:
                if not i:
                    continue
                if i.find("PL20251118")!=-1:continue  # 过滤掉这个样本
                with open(i, "r") as LEN:
                    next(LEN)  # 跳过表头
                    for len_line in LEN:
                        lines = parse_line_as_numbers(len_line)
                        if lines:
                            num_len += 1
                            text_list = lines[0:3] + [seq_type] + lines[3:]
                            df_len.loc[num_len] = text_list

            # 3) 读取 special_amplicon_stat.xls
            list_spec = os.popen(f"ls {indir}/{seq_type}/*/*/result/QC_all/*special_amplicon_stat.xls").read().strip("\n").split("\n")
            for i in list_spec:
                if not i:
                    continue
                with open(i, "r") as SPEC:
                    next(SPEC)
                    for spec_line in SPEC:
                        lines = parse_line_as_numbers(spec_line)
                        if lines:
                            num_spec += 1
                            text_list = lines[0:3] + [seq_type] + lines[3:]
                            df_spec.loc[num_spec] = text_list

            # 4) 读取 dpl_coverage_stat.xls
            list_dpl = os.popen(f"ls {indir}/{seq_type}/*/*/result/QC_all/*dpl_coverage_stat.xls").read().strip("\n").split("\n")
            for i in list_dpl:
                if not i:
                    continue
                with open(i, "r") as DPL:
                    next(DPL)
                    for dpl_line in DPL:
                        lines = parse_line_as_numbers(dpl_line)
                        if lines:
                            num_dpl += 1
                            text_list = lines[0:3] + [seq_type] + lines[3:]
                            df_dpl.loc[num_dpl] = text_list

            # 5) 读取 wes_qc_sum.txt
            list_wes = os.popen(f"ls {indir}/{seq_type}/*/*/result/QC_all/*wes_qc_sum.txt").read().strip("\n").split("\n")
            for i in list_wes:
                if not i:
                    continue
                with open(i, "r") as WES:
                    next(WES)
                    for wes_line in WES:
                        lines = parse_line_as_numbers(wes_line)
                        if lines:
                            num_wes += 1
                            df_wes.loc[num_wes] = lines

            # ---- 在这里给每个 df 添加 最小值/最大值 行 ----
            df_350= df_350.convert_dtypes()
            df_tn5= df_tn5.convert_dtypes()
            df_small= df_small.convert_dtypes()
            df_diff= df_diff.convert_dtypes()
            df_burk= df_burk.convert_dtypes()
            df_aur= df_aur.convert_dtypes()
            df_brca= df_brca.convert_dtypes()
            df_lambda= df_lambda.convert_dtypes()
            df_len= df_len.convert_dtypes()
            df_spec= df_spec.convert_dtypes()
            df_dpl= df_dpl.convert_dtypes()
            df_wes= df_wes.convert_dtypes()
            df_350 = add_min_max_rows(df_350)
            df_tn5 = add_min_max_rows(df_tn5)
            df_small = add_min_max_rows(df_small)
            df_diff = add_min_max_rows(df_diff)
            df_burk = add_min_max_rows(df_burk)
            df_aur = add_min_max_rows(df_aur)
            df_brca = add_min_max_rows(df_brca)
            df_lambda = add_min_max_rows(df_lambda)
            df_len = add_min_max_rows(df_len)
            df_spec = add_min_max_rows(df_spec)
            df_dpl = add_min_max_rows(df_dpl)
            df_wes = add_min_max_rows(df_wes)
            # ------------------------------------------

            # 将结果分别写入不同 sheet
            df_350.to_excel(writer, sheet_name='E-coli_350', index=False, header=True)
            df_tn5.to_excel(writer, sheet_name='E-coli_Tn5', index=False, header=True)
            df_small.to_excel(writer, sheet_name='E-coli_smallRNA', index=False, header=True)
            df_diff.to_excel(writer, sheet_name='E-coli_不同长度', index=False, header=True)
            df_burk.to_excel(writer, sheet_name='Burk', index=False, header=True)
            df_aur.to_excel(writer, sheet_name='金葡', index=False, header=True)
            df_brca.to_excel(writer, sheet_name='BRCA', index=False, header=True)
            df_lambda.to_excel(writer, sheet_name='lambda', index=False, header=True)
            df_len.to_excel(writer, sheet_name='length-bias', index=False, header=True)
            df_spec.to_excel(writer, sheet_name='特殊扩增子', index=False, header=True)
            df_dpl.to_excel(writer, sheet_name='DPL', index=False, header=True)
            df_wes.to_excel(writer, sheet_name='WES', index=False, header=True)

            writer.close()
    return


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-indir", required=True, help="need dir")
    parse.add_argument("-infile", required=True, help="need file")
    parse.add_argument("--outdir", default=os.getcwd(), help="output directory")
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.indir, arg.outdir)


if __name__ == "__main__":
    main()
