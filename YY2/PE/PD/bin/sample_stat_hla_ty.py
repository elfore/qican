import argparse
import os
import sys
from collections import Counter

def stat(infile, outfile, outdir,infile2):
    dicinfo = {}
    dic_t={}
    dic_n={}
    dic_A={}
    dic_B={}
    dic_C={}
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]
            dicinfo[lines[2]] = text
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
    with open(f"{outdir}/{outfile}", "w") as RAW:
        for filedir in dicinfo:
            file_list = os.popen(f"ls {filedir}/PD/result/hla_analysis_new2/*hla_result.tsv").read().strip("\n").split("\n")
            for file1 in file_list:
                pici=dicinfo[filedir].split("|")[0]
                sample_name = file1.split("/")[-1].split("_hla")[0].split("-")[0]
                sample_name2 = file1.split("/")[-1].split("_hla")[0]
                if sample_name2.find("_") != -1:
                    if pici=="PL2402201":
                        sample_name2 = sample_name2.split("_")[1]
                        sample_name = sample_name2.split("-")[0]
                    else:
                        sample_name2 = "-".join(sample_name2.split("_")[1].split("-")[1:])
                        sample_name = sample_name2.split("-")[0]
                if (sample_name2 != "CN001799-S08-D01-L09-C220C159" and sample_name2 != "CN001799-S08-D01-L10-C220C163"):
                    dic_n.setdefault(sample_name, []).append(pici + "_" + sample_name2)
                    with open(file1, "r") as FILE:
                        for line in FILE:
                            if line.startswith("Gene"):continue
                            lines = line.strip("\n").split("\t")
                            type=lines[0].split("_")[0][-1]
                            fre=":".join(lines[1].split(":")[0:2])[1:]
                            dic_t.setdefault(sample_name, {}).setdefault(type, {}).setdefault(pici + "_" + sample_name2,[]).append(fre)
        for sample_name in dic_t:
            alls=dic_n[sample_name]
            # tmpa=0
            # for sp in alls:
            #     if sp.find("SKII13126")!=-1 or sp.find("SKII13126")!=-1:
            #         tmpa+=1
            # if tmpa==0:continue
            head = sample_name +"\t"+ "一致率"+"\t" + "\t".join(alls) + "\n"
            RAW.write(head)
            for vars in dic_t[sample_name]:
                listv = []
                for name in alls:
                    if dic_t[sample_name][vars].get(name):
                        listt=dic_t[sample_name][vars][name]
                        tll=[]
                        if vars=="A":
                            for tt in listt:
                                if dic_A.get(tt):
                                    tll.append(dic_A[tt])
                        elif vars=="B":
                            for tt in listt:
                                if dic_B.get(tt):
                                    tll.append(dic_B[tt])
                        elif vars=="C":
                            for tt in listt:
                                if dic_C.get(tt):
                                    tll.append(dic_C[tt])
                        if tll==[]:
                            listv.append("阴性")
                        else:
                            tll=sorted(list(set(tll)))
                            listv.append("----".join(tll))
                    else:
                        listv.append("NA")
                count = Counter(listv)
                count_values = list(count.values())
                nums="%.2f"%((max(count_values)/len(listv))*100)+"%"
                text = "HLA-"+vars +"\t"+nums+ "\t" + "\t".join(listv) + "\n"
                RAW.write(text)
            RAW.write("\n")
    return

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-infile2", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.outfile, arg.outdir,arg.infile2)


if __name__ == "__main__":
    main()