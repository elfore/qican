import argparse
import os
import sys
import numpy as np

def stat(infile, outfile, outdir, pos):
    dicinfo = {}
    dic={}
    list_pos=[]
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]
            dicinfo[lines[2]] = text
    head1=["runID","mechineID","sampleID"]
    head2=[]
    with open(f"{pos}", "r") as pos:
        next(pos)
        for line in pos:
            lines=line.strip("\n").split("\t")
            key=lines[0]+"|"+lines[1]+"|"+lines[2]+"|"+lines[3]
            dic.setdefault(key,lines[3])
            list_pos.append(key)
            head2.append("-".join(lines[0:4]))
    head1.extend(head2)
    with open(f"{outdir}/{outfile}", "w") as RAW:
        head_text="\t".join(head1)+"\n"
        RAW.write(head_text)
        for filedir in dicinfo:
            file_list = os.popen(f"ls {filedir}/PB/result/combine_mut/*rs.vcf").read().strip("\n").split("\n")
            for file1 in file_list:
                dic_t = {}
                pici=dicinfo[filedir].split("|")[0]
                if pici=="PL2402201":
                    sample_name2=file1.split("/")[-1].split("_")[1].split(".")[0]
                else:
                    sample_name2="-".join(file1.split("/")[-1].split("_")[1].split(".")[0].split("-")[1:])
                with open(file1, "r") as FILE:
                    next(FILE)
                    for line in FILE:
                        lines = line.strip("\n").split("\t")
                        chr = lines[0]
                        poss = lines[1]
                        ref = lines[2]
                        alt = lines[3]
                        fre = lines[6]
                        if lines[6]!="-":
                            key = chr + "|" + poss + "|" + ref + "|" + alt
                            if dic.get(key) and float(fre.split("%")[0]) >= 1:
                                dic_t.setdefault(key, fre)
                m_id=pici.split("_")[1] if pici!="PL2402201" else "-"
                text_list1=[pici,m_id,sample_name2]
                text_list2=[]
                for vars in list_pos:
                    if dic_t.get(vars):
                        text_list2.append(dic_t[vars])
                    else:
                        text_list2.append("0%")
                text_list1.extend(text_list2)
                RAW.write("\t".join(text_list1)+"\n")
    return



def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-pos", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.outfile, arg.outdir,arg.pos)


if __name__ == "__main__":
    main()