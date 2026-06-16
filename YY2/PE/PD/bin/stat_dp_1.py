import argparse
import os
import sys
from datetime import datetime



def stat(infile,pos, outfile, outdir):
    dic={}
    dicinfo = {}
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]
            dicinfo[lines[2]] = text
    list_pos=[]
    head1=["runID","mechineID","sampleID"]
    head2=[]
    with open(f"{pos}", "r") as pos:
        next(pos)
        for line in pos:
            lines=line.strip("\n").split("\t")
            key=lines[0]+"|"+lines[1]
            list_pos.append(key)
            head2.append("-".join(lines[0:4]))
    head1.extend(head2)
    with open(f"{outdir}/{outfile}", "w") as RAW:
        head_text="\t".join(head1)+"\n"
        RAW.write(head_text)
        for filedir in dicinfo:
            dic_c={}
            pici=dicinfo[filedir].split("|")[0]
            cb_list=os.popen(f"ls {filedir}/PD/result2/QC/*json").read().strip("\n").split("\n")
            dp_list = os.popen(f"ls {filedir}/PD/result2/bam/*base_depth.txt").read().strip("\n").split("\n")
            for fastp in cb_list:
                sample_name2="-".join(fastp.split("/")[-1].split("_")[1].split(".")[0].split("-")[1:])
                lines = os.popen(f"head -45 {fastp}").read().strip("\n").split("\n")
                c_base=lines[17].split(":")[1].split(",")[0]
                dic_c.setdefault(sample_name2,c_base)
            for dp_file in dp_list:
                dic_dp={}
                sample_name ="-".join(dp_file.split("/")[-1].split("_")[1].split(".")[0].split("-")[1:])
                c_s=float(dic_c[sample_name])
                with open(f"{dp_file}","r") as VAR:
                    for line in VAR:
                        lines = line.strip("\n").split("\t")
                        key=lines[0]+"|"+lines[1]
                        dic_dp[key]=lines
                f_r_list=[]
                for dp_pos in list_pos:
                    if dic_dp.get(dp_pos):
                        nums=(float(dic_dp[dp_pos][2])/c_s)*10000
                    else:
                        nums=0
                    f_r_list.append(nums)
                ms=pici.split("_")[1]
                text=pici+"\t"+ms+"\t"+sample_name+"\t"+"\t".join(map(str,f_r_list))+"\n"
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