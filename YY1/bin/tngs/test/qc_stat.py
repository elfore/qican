import argparse
import glob
import os
import sys


def warn(message):
    print(f"[WARN] {message}", file=sys.stderr)


def nonempty_files(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        warn(f"no files matched: {pattern}")
    return files


def has_header(path):
    if os.path.getsize(path) == 0:
        warn(f"skip empty file: {path}")
        return False
    return True


def stat(infile,pos, outdir):
    dicinfo = {}
    dic_pos={}
    dic={}
    dic_dp={}
    dic_qc={}
    with open(f"{pos}", "r") as pos:
        for line in pos:
            lines=line.strip("\n").split("\t")
            rl=lines[1].split(",")
            dic_pos.setdefault(lines[0],rl)
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[1] + "|" + lines[2]
            dicinfo[lines[3]] = text
    for filedir in dicinfo:
        pici=dicinfo[filedir].split("|")[0]
        file_list = nonempty_files(f"{filedir}/tNGS/LC*/LC*raw_result.xls")
        file_list2 = nonempty_files(f"{filedir}/tNGS/LC*/LC*qc.xls")
        for file1 in file_list:
            if not has_header(file1):
                continue
            sample_name=file1.split("/")[-1].split(".")[0]
            if sample_name not in dic_pos:
                warn(f"skip sample without result_info entry: {sample_name}")
                continue
            with open(file1, "r") as FILE:
                try:
                    next(FILE)
                except StopIteration:
                    warn(f"skip headerless file: {file1}")
                    continue
                list_tp=[]
                for line in FILE:
                    lines = line.strip("\n").split("\t")
                    ty=lines[1]
                    dic_dp.setdefault(sample_name, {}).setdefault(pici, {}).setdefault(lines[1], lines[6]+"("+lines[4]+")")
                    if lines[1]=="vanA" or lines[1]=="热带念珠菌":continue
                    if lines[6].find("Pass") != -1:
                        list_tp.append(lines[1])
                        if ty in dic_pos[sample_name]:
                            dic.setdefault(pici,{}).setdefault(sample_name, {}).setdefault("zy", []).append(lines[1])
                        else:
                            dic.setdefault(pici,{}).setdefault(sample_name, {}).setdefault("jy", []).append(lines[1])
                for res in dic_pos[sample_name]:
                    if res not in list_tp:
                        dic.setdefault(pici,{}).setdefault(sample_name, {}).setdefault("lj", []).append(res)
        for file2 in file_list2:
            if not has_header(file2):
                continue
            sample_name = file2.split("/")[-1].split(".")[0]
            with open(file2, "r") as FILE:
                try:
                    next(FILE)
                except StopIteration:
                    warn(f"skip headerless file: {file2}")
                    continue
                for line in FILE:
                    dic_qc.setdefault(pici, {}).setdefault(sample_name, line.strip("\n"))
    with open(f"{outdir}/tngs_stat.txt", "w") as RAW:
        headers = [
            "批次",
            "仪器id",
    "sampleID",
    "raw_reads",
    "clean_reads",
    "adapter(%)",
    "raw_length",
    "clean_length",
    "raw_Q20(%)",
    "mapping_reads",
    "mapping_rate(%)",
    "uniq_mapping_rate(%)",
    "multi_mapping_reads",
    "multi_mapping_rate(%)",
    "amplicon>1",
    "amplicon>20",
    "amplicon>100",
    "amplicon>500",
    "真阳性结果",
    "假阳性结果",
    "漏检结果",
]
        RAW.write("\t".join(headers) + "\n")
        for pici in dic:
            for sample_name in dic[pici]:
                yq_id = pici.split("_")[1]
                text_list = [pici, yq_id]
                qc_line = dic_qc.get(pici, {}).get(sample_name)
                if qc_line is None:
                    warn(f"skip sample without qc entry: {pici} {sample_name}")
                    continue
                qc_info = qc_line.strip("\n").split("\t")
                text_list.extend(qc_info)
                zy = dic[pici][sample_name].get("zy", [])
                jy = dic[pici][sample_name].get("jy", [])
                lj = dic[pici][sample_name].get("lj", [])
                zy_list=[zy_item + "_" + dic_dp[sample_name][pici][zy_item] if zy_item in dic_dp[sample_name][pici] else zy_item for zy_item in zy]
                jy_list=[jy_item + "_" + dic_dp[sample_name][pici][jy_item] if jy_item in dic_dp[sample_name][pici] else jy_item for jy_item in jy]
                lj_list=[lj_item + "_" + dic_dp[sample_name][pici][lj_item] if lj_item in dic_dp[sample_name][pici] else lj_item for lj_item in lj]
                text_list.append(",".join(zy_list))
                text_list.append(",".join(jy_list))
                text_list.append(",".join(lj_list))
                text_list = [str(item) for item in text_list]
                RAW.write("\t".join(text_list) + "\n")
    return



def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-pos", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile,  arg.pos,arg.outdir)


if __name__ == "__main__":
    main()
