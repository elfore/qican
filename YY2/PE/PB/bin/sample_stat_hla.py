import argparse
import os
import sys

def stat(infile, outfile, outdir):
    dicinfo = {}
    dic_t={}
    dic_n={}
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[0] + "|" + lines[1]
            dicinfo[lines[2]] = text
    with open(f"{outdir}/{outfile}", "w") as RAW:
        for filedir in dicinfo:
            file_list = os.popen(f"ls {filedir}/PB/result/hla_analysis_new2/*hla_result.tsv").read().strip("\n").split("\n")
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
                dic_n.setdefault(sample_name, []).append(pici + "_" + sample_name2)
                with open(file1, "r") as FILE:
                    for line in FILE:
                        if line.startswith("Gene"):continue
                        lines = line.strip("\n").split("\t")
                        type=lines[0].split("_")[0][-1]+lines[0][-1]
                        fre=lines[1]
                        dic_t.setdefault(sample_name, {}).setdefault(type, {}).setdefault(pici + "_" + sample_name2,fre)
        for sample_name in dic_t:
            alls=dic_n[sample_name]
            head = sample_name + "\t" + "\t".join(alls) + "\n"
            RAW.write(head)
            for vars in dic_t[sample_name]:
                listv = []
                for name in alls:
                    if dic_t[sample_name][vars].get(name):
                        listv.append(dic_t[sample_name][vars][name])
                    else:
                        listv.append("NA")
                text = vars + "\t" + "\t".join(listv) + "\n"
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