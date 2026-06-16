import argparse
import os
import sys
import numpy as np

def stat(infile, outfile, outdir):
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
            file1_list=os.popen(f"ls {filedir}/PB/result/qc_merge/*qc_merge.csv").read().strip("\n").split("\n")
            for spt in file1_list:
                with open(spt,"r") as QC:
                    next(QC)
                    for line in QC:
                        lines=line.strip("\n").split(",")
                        if pici=="PL2402201":
                            key=lines[0].split("_")[1]
                        else:
                            key="-".join(lines[0].split("-")[1:])
                        dic_s.setdefault(key,[]).extend(lines[1:])
            for key in dic_s:
                if pici.find("PL2402201")!=-1:
                    ms="-"
                else:
                    ms=pici.split("_")[1]
                text=pici+"\t"+ms+"\t"+key+"\t"+"\t".join(map(str,dic_s[key]))+"\n"
                RAW.write(text)
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