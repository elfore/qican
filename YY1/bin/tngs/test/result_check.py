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


def stat(infile, outfile, outdir, pos):
    dicinfo = {}
    dic={}
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[1] + "|" + lines[2]
            dicinfo[lines[3]] = text
    head1=["runID","mechineID","真阳性结果","假阳性结果"]
    with open(f"{pos}", "r") as pos:
        for line in pos:
            lines=line.strip("\n").split("\t")
            rl=lines[1].split(",")
            dic.setdefault(lines[0],rl)
    with open(f"{outdir}/{outfile}", "w") as RAW:
        head_text="\t".join(head1)+"\n"
        RAW.write(head_text)
        for filedir in dicinfo:
            pici=dicinfo[filedir].split("|")[0]
            m_id=pici.split("_")[1]
            file_list = nonempty_files(f"{filedir}/tNGS/LC*/LC*raw_result.xls")
            all_p=0
            true_p=0
            a_true_p=0
            for file1 in file_list:
                if not has_header(file1):
                    continue
                sample_name=file1.split("/")[-1].split(".")[0]
                if sample_name not in dic:
                    warn(f"skip sample without result_info entry: {sample_name}")
                    continue
                s_re=dic[sample_name]
                a_true_p+=len(s_re)
                with open(file1, "r") as FILE:
                    try:
                        next(FILE)
                    except StopIteration:
                        warn(f"skip headerless file: {file1}")
                        continue
                    p_list=[]
                    for line in FILE:
                        lines = line.strip("\n").split("\t")
                        if lines[1]=="vanA" or lines[1]=="热带念珠菌":continue
                        if lines[6].find("Pass")!=-1:
                            all_p+=1
                            if lines[1] in s_re:
                                true_p+=1
                                p_list.append(lines[1])
            SEN=true_p/a_true_p if a_true_p else 0
            PPV=true_p/all_p if all_p else 0
            if not a_true_p or not all_p:
                warn(f"{pici} has zero denominator: expected={a_true_p}, pass={all_p}")
            text_list=[pici,m_id,SEN,PPV]
            RAW.write("\t".join(map(str,text_list))+"\n")
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
