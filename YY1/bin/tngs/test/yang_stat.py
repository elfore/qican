import argparse
import os
import sys

def stat(infile, outdir):
    dicinfo = {}
    dic={}
    dic2={}
    dic_dp={}
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[1] + "|" + lines[2]
            dicinfo[lines[3]] = text
    for filedir in dicinfo:
        pici=dicinfo[filedir].split("|")[0]
        file_list = os.popen(f"ls {filedir}/tNGS/LC*/LC*raw_result.xls").read().strip("\n").split("\n")
        for file1 in file_list:
            sample_name=file1.split("/")[-1].split(".")[0]
            with open(file1, "r") as FILE:
                next(FILE)
                for line in FILE:
                    lines = line.strip("\n").split("\t")
                    dic_dp.setdefault(sample_name, {}).setdefault(pici, {}).setdefault(lines[1], lines[4])
                    if lines[1] == "vanA":
                        continue
                    if lines[6].find("Pass") != -1 or lines[6].find("灰区") != -1:
                        dic.setdefault(sample_name,{}).setdefault(pici, {}).setdefault(lines[1], lines[6])
                        dic2.setdefault(sample_name, []).append(lines[1])
    for sample_name in dic:
        with open(f"{outdir}/{sample_name}.txt", "w") as RAW:
            list1=list(set(dic2[sample_name]))
            head1="#批次"+"\t"+"样品名\t"+"\t".join(list1)+"\n"
            RAW.write(head1)
            for pici in dic[sample_name]:
                text_list = [pici, sample_name]
                for item in list1:
                    dp=dic_dp[sample_name][pici][item] if item in dic_dp[sample_name][pici] else "0"
                    if item in dic[sample_name][pici]:
                        text_list.append(dic[sample_name][pici][item]+"("+dp+")")
                    else:
                        text_list.append("未检出"+"("+dp+")")
                RAW.write("\t".join(text_list) + "\n")
    return



def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile,arg.outdir)


if __name__ == "__main__":
    main()