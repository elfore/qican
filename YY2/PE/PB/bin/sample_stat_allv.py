import argparse
import os
import sys
import numpy as np


def get_ty(ty):
    if ty=="0|0":
        return "阴性"
    else:
        return "阳性"


def stat(infile,infile2, outfile, outdir, pos):
    dicinfo = {}
    dic_r={}
    dic_r2={}
    dic_r3={}
    dic={}
    rz_dir=""
    dic_n={}
    list_pos=[]
    dic_A={}
    dic_B={}
    dic_C={}
    with open(f"{infile2}", "r") as IN2:
        for line in IN2:
            lines = line.strip("\n").split("\t")
            text=lines[1]+"|"+lines[5]
            if lines[3]=="HLA-A":
                dic_A.setdefault(lines[0],text)
            elif lines[3]=="HLA-B":
                dic_B.setdefault(lines[0], text)
            else:
                dic_C.setdefault(lines[0], text)
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]
            if lines[0]=="PL2402201":
                rz_dir+=lines[2]
            else:
                dicinfo[lines[2]] = text
    with open(f"{pos}", "r") as pos:
        next(pos)
        for line in pos:
            lines=line.strip("\n").split("\t")
            key=lines[0]+"|"+lines[1]+"|"+lines[2]+"|"+lines[3]
            dic.setdefault(key,lines[3])
            list_pos.append(key)
    head1=["runID","mechineID","sampleID","阳性符合率","阴性符合率"]
    with open(f"{outdir}/{outfile}", "w") as RAW:
        RAW.write("\t".join(head1)+"\n")
        dic_r_num={}
        file_list = os.popen(f"ls {rz_dir}/PB/result/combine_mut/*rs.vcf").read().strip("\n").split("\n")
        file_list2 = os.popen(f"ls {rz_dir}/PB/result/hla_analysis_new2/*hla_result.tsv").read().strip("\n").split("\n")
        for file1 in file_list:
            sample_name2=file1.split("/")[-1].split("_")[1].split(".")[0]
            r_yang=0
            r_ying=0
            with open(file1, "r") as FILE:
                for line in FILE:
                    if line.startswith("chrom"):continue
                    lines = line.strip("\n").split("\t")
                    chr = lines[0]
                    poss = lines[1]
                    ref = lines[2]
                    alt = lines[3]
                    ty = lines[7]
                    key = chr + "|" + poss + "|" + ref + "|" + alt
                    if dic.get(key):
                        ty_fre=get_ty(ty)
                        if ty_fre=="阴性":
                            r_ying+=1
                        else:
                            r_yang+=1
                        dic_r.setdefault(sample_name2, {}).setdefault(key, ty_fre)
                for key in list_pos:
                    if not dic_r[sample_name2].get(key):
                        dic_r.setdefault(sample_name2, {}).setdefault(key, "阴性")
                        r_ying+=1
            dic_r_num.setdefault(sample_name2, {}).setdefault("ying", r_ying)
            dic_r_num.setdefault(sample_name2, {}).setdefault("yang", r_yang)
        for file2 in file_list2:
            sample_name2 =file2.split("/")[-1].split("_hla")[0].split("_")[1]
            with open(file2, "r") as FILE:
                tll_dic={"HLA-A":[],"HLA-B":[],"HLA-C":[]}
                yyy=""
                for line in FILE:
                    if line.startswith("Gene"):continue
                    lines = line.strip("\n").split("\t")
                    type=lines[0].split("_")[0]
                    fre=lines[1]
                    fre2=":".join(lines[1].split(":")[0:2])[1:]
                    if type=="HLA-A":
                        if dic_A.get(fre2):
                            tll_dic.setdefault("HLA-A",[]).append(dic_A[fre2])
                    elif type=="HLA-B":
                        if dic_B.get(fre2):
                            tll_dic.setdefault("HLA-B",[]).append(dic_B[fre2])
                    elif type=="HLA-C":
                        if dic_C.get(fre2):
                            tll_dic.setdefault("HLA-C",[]).append(dic_C[fre2])
                    dic_r2.setdefault(sample_name2, {}).setdefault(type, []).append(fre)
                for tyy in tll_dic:
                    tll=tll_dic[tyy]
                    if tll==[]:
                        yyy="阴性"
                    else:
                        tll=sorted(list(set(tll)))
                        yyy=("----".join(tll))
                    dic_r3.setdefault(sample_name2, {}).setdefault(tyy, yyy)
        for filedir in dicinfo:
            file_list = os.popen(f"ls {filedir}/PB/result/combine_mut/*rs.vcf").read().strip("\n").split("\n")
            file_list2 = os.popen(f"ls {filedir}/PB/result/hla_analysis_new2/*hla_result.tsv").read().strip("\n").split("\n")
            dic_f={}
            for file1 in file_list:
                pici=dicinfo[filedir].split("|")[0]
                sample_name = file1.split("/")[-1].split("_")[1].split(".")[0].split("-")[1]
                sample_name2="-".join(file1.split("/")[-1].split("_")[1].split(".")[0].split("-")[1:])
                ying=0
                yang=0
                dic_t = {}
                dic_n.setdefault(sample_name,[]).append(pici+"_"+sample_name2)
                with open(file1, "r") as FILE:
                    for line in FILE:
                        if line.startswith("chrom"):continue
                        lines = line.strip("\n").split("\t")
                        chr = lines[0]
                        poss = lines[1]
                        ref = lines[2]
                        alt = lines[3]
                        ty = lines[7]
                        ty_fre=get_ty(ty)
                        key = chr + "|" + poss + "|" + ref + "|" + alt
                        if dic.get(key):
                            dic_t.setdefault(key, ty_fre)
                            if ty_fre=="阳性" and ty_fre==dic_r[sample_name2][key]:
                                yang+=1
                            if ty_fre=="阴性" and ty_fre==dic_r[sample_name2][key]:
                                ying+=1
                    for key in list_pos:
                        if not dic_t.get(key):
                            if dic_r[sample_name2][key]=="阴性":
                                ying+=1
                yang_result=yang/dic_r_num[sample_name2]["yang"]
                ying_result=ying/dic_r_num[sample_name2]["ying"]
                m_id=pici.split("_")[1]
                text_list=[pici,m_id,sample_name2,yang_result,ying_result]
                dic_f.setdefault(sample_name2, []).extend(text_list)
            for file2 in file_list2:
                sample_name2 ="-".join(file2.split("/")[-1].split("_hla")[0].split("_")[1].split("-")[1:])
                c_num=0
                dic_tt={}
                yy_result=""
                tll_dic={"HLA-A":[],"HLA-B":[],"HLA-C":[]}
                cc_num=0
                with open(file2, "r") as FILE:
                    for line in FILE:
                        if line.startswith("Gene"):continue
                        lines = line.strip("\n").split("\t")
                        type=lines[0].split("_")[0]
                        fre=lines[1]
                        fre2=":".join(lines[1].split(":")[0:2])[1:]
                        if type=="HLA-A":
                            if dic_A.get(fre2):
                                tll_dic.setdefault("HLA-A",[]).append(dic_A[fre2])
                        elif type=="HLA-B":
                            if dic_B.get(fre2):
                                tll_dic.setdefault("HLA-B",[]).append(dic_B[fre2])
                        elif type=="HLA-C":
                            if dic_C.get(fre2):
                                tll_dic.setdefault("HLA-C",[]).append(dic_C[fre2])
                        dic_tt.setdefault(type, []).append(fre)
                for tyy in tll_dic:
                    tll=tll_dic[tyy]
                    if tll==[]:
                        yy_result="阴性"
                    else:
                        tll=sorted(list(set(tll)))
                        yy_result=("----".join(tll))
                    true_result=dic_r3[sample_name2][tyy]
                    if yy_result==true_result:
                        cc_num+=1
                for ty in dic_tt:
                    ty_list=dic_tt[ty]
                    if not ty_list[0] in dic_r2[sample_name2][ty]:
                        c_num+=1
                    if not ty_list[1] in dic_r2[sample_name2][ty]:
                        c_num+=1
                    if ty_list[0]==ty_list[1] and (ty_list[0] in dic_r2[sample_name2][ty]) and (dic_r2[sample_name2][ty][0]!=dic_r2[sample_name2][ty][1]):
                        c_num+=1
                
                hla_result=(6-c_num)/6
                hla_result2=cc_num/3
                r1=dic_f[sample_name2]
                r1.append(hla_result)
                r1.append(hla_result2)
                RAW.write("\t".join(map(str,r1))+"\n")
    return


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-infile2", required=True, help="need dir")
    parse.add_argument("-pos", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile,arg.infile2, arg.outfile, arg.outdir,arg.pos)


if __name__ == "__main__":
    main()