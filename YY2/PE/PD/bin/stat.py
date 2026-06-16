import argparse
import os
import sys
import numpy as np


# "header": runID	mechineID	sampleID	raw_reads	clean_reads	effective_rate	adapter	raw_len	filter_len	Q20	Q30	GC	mapped_reads	mapped_rate	coverage	avg_depth	on_target_whole	on_target_hotspot	on_target_cyp	on_target_hla_a	on_target_hla_b	uniformity_snp	uniformity_lg	hotspot_cov	hotspot_50_cov	pb_cov	pb_50_cov	hla_a_cov	hla_a_50_cov	hla_b_cov	hla_b_50_cov	hla_c_cov	hla_c_50_cov	var_num	un_p_r	cov	cov_20	cov_30	cov_50	cov_100	cov_200	cov_500	cov_1000	是否合格
def get_pos_stat(infile,list,file_list2):
    dic={}
    f_list=[]
    dic_r={}
    with open(f"{infile}", "r") as pos:
        for line in pos:
            lines=line.strip("\n").split("\t")
            text=lines[3]
            dic.setdefault(lines[0],{}).setdefault(lines[1],[]).append(text)
    for qcr in list:
        file = qcr.split("/")[-1]
        sample_name =file.split(".")[0]
        a=0;b=0;c=0;d=0;e=0;f=0;all=0;ff=0;bb=0;bbb=0
        with open(f"{qcr}","r") as VAR:
            for line in VAR:
                lines = line.strip("\n").split("\t")
                if dic.get(lines[0]):
                    if dic[lines[0]].get(lines[1]):
                        all+=1
                        if int(lines[2])>0:
                            a+=1
                        else:
                            ff+=len(dic[lines[0]][lines[1]])
                        if int(lines[2])>=20:
                            b+=1
                        if int(lines[2])>=30:
                            bbb+=1
                        if int(lines[2])>=50:
                            bb+=1
                        if int(lines[2])>=100:
                            c+=1
                        if int(lines[2])>=200:
                            d+=1
                        if int(lines[2])>=500:
                            e+=1
                        if int(lines[2])>=1000:
                            f += 1
        if (bbb/all)==1:
            zk="合格"
        else:
            zk="不合格"
        text=sample_name+"\t"+str(all)+"\t"+str(ff)+"\t"+str(round((a/all)*100,2))+"%("+str(a)+")\t"+str(round((b/all)*100,2))+"%("+str(b)+")\t"+str(round((bbb/all)*100,2))+"%("+str(bbb)+")\t"+str(round((bb/all)*100,2))+"%("+str(bb)+")\t"+str(round((c/all)*100,2))+"%("+str(c)+")\t"+str(round((d/all)*100,2))+"%("+str(d)+")\t"+str(round((e/all)*100,2))+"%("+str(e)+")\t"+str(round((f/all)*100,2))+"%("+str(f)+")"+"\t"+zk
        dic_r.setdefault(sample_name,text)
    for file2 in file_list2:
        sample_name2 =file2.split("/")[-1].split("_hla")[0]
        a=0
        with open(file2, "r") as FILE:
            for line in FILE:
                if line.startswith("Gene"):continue
                lines = line.strip("\n").split("\t")
                if float(lines[2])<0.75:
                    print(a)
                    a+=1
        b=dic_r[sample_name2].split("\t")[-1]
        if b=="合格" and a==0:
            f_list.append(dic_r[sample_name2]+"\n")
        else:
            f_list.append("\t".join(dic_r[sample_name2].split("\t")[:-1])+"\t"+"不合格"+"\n")

    return f_list



def stat(infile,pos, outfile, outdir):
    dicinfo = {}
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]
            dicinfo[lines[2]] = text
    with open(f"{outdir}/{outfile}", "w") as RAW:
        for filedir in dicinfo:
            pici=dicinfo[filedir].split("|")[0]
            dic_s={}
            file1=os.popen(f"ls {filedir}/PD/result/QC/*final_QC_stat.xls").read().strip("\n").split("\n")[0]
            file2=os.popen(f"ls {filedir}/PD/result/bam_stat/*final_bam_stat.xls").read().strip("\n").split("\n")[0]
            dp_list = os.popen(f"ls {filedir}/PD/result/base_depth/*base_depth.txt").read().strip("\n").split("\n")
            hla_list=os.popen(f"ls {filedir}/PD/result/hla_analysis_new2/*hla_result.tsv").read().strip("\n").split("\n")
            pos_dp_qc=get_pos_stat(pos,dp_list,hla_list)
            with open(file1,"r") as QC,open(file2,"r") as BAM:
                next(QC)
                next(BAM)
                for line in QC:
                    lines=line.strip("\n").split("\t")
                    if pici=="PL2402201":
                        key=lines[0].split("_")[1]
                    else:
                        key="-".join(lines[0].split("-")[1:])
                    dic_s.setdefault(key,[]).extend(lines[1:])
                for line in BAM:
                    lines=line.strip("\n").split("\t")
                    if pici=="PL2402201":
                        key=lines[0].split("_")[1]
                    else:
                        key="-".join(lines[0].split("-")[1:])
                    dic_s.setdefault(key,[]).extend(lines[1:])
                for line in pos_dp_qc:
                    lines=line.strip("\n").split("\t")
                    if pici=="PL2402201":
                        key=lines[0].split("_")[1]
                    else:
                        key="-".join(lines[0].split("-")[1:])
                    dic_s.setdefault(key,[]).extend(lines[1:])
            for key in dic_s:
                if pici=="PL2402201":
                    ms="-"
                else:
                    ms=pici.split("_")[1]
                text=pici+"\t"+ms+"\t"+key+"\t"+"\t".join(map(str,dic_s[key]))+"\n"
                RAW.write(text)
    return


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-pos", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.pos,arg.outfile, arg.outdir)


if __name__ == "__main__":
    main()