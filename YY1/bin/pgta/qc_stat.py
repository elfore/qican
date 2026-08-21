import argparse
import os
import sys
import re




def extract_chromosome_and_p_positions(input_string):
    # 正则表达式分别匹配染色体（支持 chr 或 q 型）
    pattern_chr = r'([+-]?chr[\d]+|[+-]?q[\d]+)'  # 支持数字和字母的染色体类型（chr 和 q）
    match_chr = re.search(pattern_chr, input_string)
    
    # 正则表达式提取 p 位置部分
    pattern_p = r'(p\d+\.\d+|q\d+\.\d+)'  # 匹配 p 或 q 型的位置
    match_p = re.findall(pattern_p, input_string)
    
    if match_chr and match_p:
        # 提取染色体部分并去掉 + 或 - 符号
        chromosome = match_chr.group(1)
        chromosome = re.sub(r'^[+-]', '', chromosome)  # 去掉开头的 + 或 -
        
        # 提取第一个 p 位置
        p_position_1 = match_p[0]
        # 提取第二个 p 位置（如果存在）
        p_position_2 = match_p[1] if len(match_p) > 1 else None
        
        return chromosome, p_position_1, p_position_2
    else:
        return None, None, None


def check_vars(vars,t_result,tp):
    if tp=="假阳":
        f_t="阴性"
    else:
        f_t="未见异常"
    if t_result.find(f_t)!=-1:
        vars=f"{tp}："+vars
    else:
        if vars.find("chrX")!=-1 or vars.find("chrY")!=-1:
            return vars
        else:
            a=0
            t_list=t_result.split(";")
            for t in t_list:
                t_chr,t_p1,t_p2=extract_chromosome_and_p_positions(t)
                tt=[t_p1,t_p2]
                v_chr,v_p1,v_p2=extract_chromosome_and_p_positions(vars)
                if v_chr==t_chr and (v_p1 in tt or v_p2 in tt):
                    a+=1
            if a==0:
                vars=f"{tp}："+vars
                return vars
    return vars

def check_fields_in_sublists(all_list):
    # 定义需要检查的字段
    target_fields = ['漏检', '假阳']
    
    # 用于存储找到的字段
    found_fields = set()
    # 遍历大列表中的每个子列表
    for sublist in all_list:
        for item in sublist:
            # 如果字段在当前项中，加入到 found_fields 集合中
            for field in target_fields:
                if field in item:
                    found_fields.add(field)
    
    # 根据找到的字段返回相应的结果
    if '漏检' in found_fields and '假阳' in found_fields:
        return '漏检, 假阳'
    elif '漏检' in found_fields:
        return '漏检'
    elif '假阳' in found_fields:
        return '假阳'
    else:
        return '一致'

def stat(infile, infile2,outfile, outdir):
    dicinfo = {}
    diccnv={}
    with open(f"{infile}", "r") as IN,open(infile2,"r") as IN2:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]+ "|" + lines[2]
            dicinfo[lines[3]] = text
        for line in IN2:
            lines = line.strip("\n").split("\t")
            diccnv[lines[0]]=lines[1]
    with open(f"{outdir}/{outfile}", "w") as RAW:
        for filedir in dicinfo:
            dic={}
            pici=dicinfo[filedir].split("|")[1]
            yiqi=dicinfo[filedir].split("|")[1].split("_")[1]
            tp=dicinfo[filedir].split("|")[2]
            json_list=os.popen(f"ls {filedir}/PGT/*/*json").read().strip("\n").split("\n")
            for fastp in json_list:
                sample_name = fastp.split("/")[-2]
                lines = os.popen(f"head -45 {fastp}").read().strip("\n").split("\n")
                r_20=lines[9].split(":")[1].split(",")[0]
                r_30=lines[10].split(":")[1].split(",")[0]
                dic.setdefault(sample_name, {}).setdefault('r_20', r_20)
                dic.setdefault(sample_name, {}).setdefault('r_30', r_30)
            file_list = os.popen(f"ls {filedir}/PGT/*/*totalqc.txt").read().strip("\n").split("\n")
            for file1 in file_list:
                sample_name = file1.split("/")[-2]
                with open(file1, "r") as FILE:
                    for line in FILE:
                        if line.startswith("sampleID"):continue
                        lines=line.strip("\n").split("\t")
                        r_reads=lines[1]
                        u_reads=r_reads if float(r_reads)<2000000 else 2000000
                        m_rate=lines[18]
                        u_m_rate=lines[19]
                        m_q30_rate=lines[20]
                        mb_cv=lines[17]
                        g_cov=lines[21]
                        gd=lines[14]
                        result=lines[22]
                    r_list=result.strip(";").split(";")[1:]
                    cnv_4_list1=[];cnv_4_list2=[];mos_4_list1=[];mos_4_list2=[]
                    t_result=diccnv[sample_name]
                    if result.find("未见异常")!=-1:
                        pass
                    else:         
                        for vars in r_list:
                            if vars.find("chromsome")!=-1:continue
                            if vars.find("mos")!=-1:
                                if float(vars.split(",")[1].split("M")[0])<4:
                                    vars=check_vars(vars,t_result,"假阳")
                                    mos_4_list1.append(vars)
                                else:
                                    vars=check_vars(vars,t_result,"假阳")
                                    mos_4_list2.append(vars)
                            else:
                                if float(vars.split(",")[1].split("M")[0])<4:
                                    vars=check_vars(vars,t_result,"假阳")
                                    cnv_4_list1.append(vars)
                                else:
                                    vars=check_vars(vars,t_result,"假阳")
                                    cnv_4_list2.append(vars)
                    if t_result=="阴性":
                        pass
                    else:
                        r_result=";".join(r_list)
                        for vars in t_result.split(";"):
                            if vars.find("mos")!=-1:
                                if float(vars.split(",")[1].split("M")[0])<4:
                                    vars=check_vars(vars,r_result,"漏检")
                                    if vars.find("漏检")!=-1:
                                        mos_4_list1.append(vars)
                                else:
                                    vars=check_vars(vars,r_result,"漏检")
                                    if vars.find("漏检")!=-1:
                                        mos_4_list2.append(vars)
                            else:
                                if float(vars.split(",")[1].split("M")[0])<4:
                                    vars=check_vars(vars,r_result,"漏检")
                                    if vars.find("漏检")!=-1:
                                        cnv_4_list1.append(vars)
                                else:
                                    vars=check_vars(vars,r_result,"漏检")
                                    if vars.find("漏检")!=-1:
                                        cnv_4_list2.append(vars)
                    all_list=[mos_4_list1,mos_4_list2,cnv_4_list1,cnv_4_list2]
                    f_result=check_fields_in_sublists(all_list)
                    text_list=[pici,yiqi,tp,"","",sample_name,t_result,r_reads,dic[sample_name]["r_20"],dic[sample_name]["r_30"],u_reads,m_rate,u_m_rate,m_q30_rate,mb_cv,g_cov,gd,";".join(cnv_4_list1),";".join(cnv_4_list2),";".join(mos_4_list1),";".join(mos_4_list2),f_result]
                    text="\t".join(map(str,text_list))+"\n"
                    RAW.write(text)
    return
    
def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-infile2", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.infile2,arg.outfile, arg.outdir)


if __name__ == "__main__":
    main()