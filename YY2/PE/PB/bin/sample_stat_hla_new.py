import argparse
import os
import sys

def stat(infile, outfile, outdir):
    dicinfo = {}
    dic_t={}
    dic_n={}
    
    # 1. 读取基础信息
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]
            dicinfo[lines[2]] = text
            
    # 2. 解析每个样本的鉴定结果
    with open(f"{outdir}/{outfile}", "w") as RAW:
        for filedir in dicinfo:
            file_list = os.popen(f"ls {filedir}/PB/result/hla_analysis_new2/*hla_result.tsv").read().strip("\n").split("\n")
            if not file_list or file_list[0] == "":
                continue
                
            for file1 in file_list:
                pici = dicinfo[filedir].split("|")[0]
                sample_name = file1.split("/")[-1].split("_hla")[0].split("-")[0]
                sample_name2 = file1.split("/")[-1].split("_hla")[0]
                
                if sample_name2.find("_") != -1:
                    if pici == "PL2402201":
                        sample_name2 = sample_name2.split("_")[1]
                        sample_name = sample_name2.split("-")[0]
                    else:
                        sample_name2 = "-".join(sample_name2.split("_")[1].split("-")[1:])
                        sample_name = sample_name2.split("-")[0]
                        
                test_id = pici + "_" + sample_name2
                
                # 为了防止重复添加列名，加入判断
                if sample_name not in dic_n:
                    dic_n[sample_name] =[]
                if test_id not in dic_n[sample_name]:
                    dic_n[sample_name].append(test_id)
                    
                with open(file1, "r") as FILE:
                    for line in FILE:
                        if line.startswith("Gene"): continue
                        lines = line.strip("\n").split("\t")
                        if len(lines) < 2: continue
                        type = lines[0].split("_")[0][-1] + lines[0][-1]
                        fre = lines[1]
                        dic_t.setdefault(sample_name, {}).setdefault(type, {})[test_id] = fre

        # 3. 对齐数据、计算一致率并输出
        for sample_name in dic_t:
            alls = dic_n[sample_name]
            if not alls: continue
            
            baseline = alls[0] # 以每种样品的第一列样品作为基准
            tests = alls[1:]   # 后续要计算一致率的样品
            
            # --- 新增功能 1：对齐 A1/A2, B1/B2 (解决颠倒问题) ---
            # 提取基因类型 (例如从 'A1', 'A2' 中提取 'A')
            genes = set([var[:-1] for var in dic_t[sample_name].keys()]) 
            for gene in genes:
                v1 = gene + "1" # 如 A1
                v2 = gene + "2" # 如 A2
                
                # 如果A1和A2都有，则进行匹配对齐
                if v1 in dic_t[sample_name] and v2 in dic_t[sample_name]:
                    b1 = dic_t[sample_name][v1].get(baseline, "NA")
                    b2 = dic_t[sample_name][v2].get(baseline, "NA")
                    
                    for t in tests:
                        t1 = dic_t[sample_name][v1].get(t, "NA")
                        t2 = dic_t[sample_name][v2].get(t, "NA")
                        
                        # 方案1(不交换)的分数
                        score_no_swap = (1 if t1 == b1 and b1 != "NA" else 0) + (1 if t2 == b2 and b2 != "NA" else 0)
                        # 方案2(交换后)的分数
                        score_swap = (1 if t2 == b1 and b1 != "NA" else 0) + (1 if t1 == b2 and b2 != "NA" else 0)
                        
                        # 如果交换后的匹配度更高，就在字典里颠倒它们的位置
                        if score_swap > score_no_swap:
                            dic_t[sample_name][v1][t] = t2
                            dic_t[sample_name][v2][t] = t1

            # --- 新增功能 2 & 3：计算一致率和输出排版 ---
            # 写入表头，插入 "一致率" 列
            head = sample_name + "\t" + "一致率\t" + "\t".join(alls) + "\n"
            RAW.write(head)
            
            # 使用 sorted 确保按 A1, A2, B1, B2 的顺序输出
            for vars in sorted(dic_t[sample_name].keys()):
                listv = []
                b_val = dic_t[sample_name][vars].get(baseline, "NA") # 获取基准值
                match_count = 0
                total_compare = len(tests)
                
                for name in alls:
                    val = dic_t[sample_name][vars].get(name, "NA")
                    listv.append(val)
                    
                    # 如果不是基准列，且不为空，参与匹配统计
                    if name != baseline:
                        if val == b_val and b_val != "NA":
                            match_count += 1
                            
                # 计算一致率
                if total_compare > 0:
                    rate_str = f"{int((match_count / total_compare) * 100)}%"
                else:
                    rate_str = "NA" # 如果后面没有样本可供对比
                    
                # 拼接输出行：指标名(如A1) + 一致率(如100%) + 各列分型结果
                text = vars + "\t" + rate_str + "\t" + "\t".join(listv) + "\n"
                RAW.write(text)
                
            RAW.write("\n")
    return

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.outfile, arg.outdir)

if __name__ == "__main__":
    main()