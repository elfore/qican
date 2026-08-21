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


def stat(infile, outfile, outdir, pos,pos2):
    dicinfo = {}
    list_id=[]
    dic_sp={}
    with open(f"{infile}", "r") as IN:
        for line in IN:
            lines = line.strip("\n").split("\t")
            text = lines[1] + "|" + lines[2]
            dicinfo[lines[3]] = text
    with open(f"{pos}", "r") as pos:
        next(pos)
        for line in pos:
            id=line.strip("\n")
            list_id.append(id)
    with open(f"{pos2}", "r") as pos2:
        for line in pos2:
            lines = line.strip("\n").split("\t")
            dic_sp[lines[0]]=lines[1]+"|"+lines[2]
    head1=["runID","mechineID","sampleID"]
    head1.extend(list_id)
    with open(f"{outdir}/{outfile}", "w") as RAW:
        head_text="\t".join(head1)+"\n"
        RAW.write(head_text)
        for filedir in dicinfo:
            fz=0
            fm=0
            pici=dicinfo[filedir].split("|")[0]
            m_id=pici.split("_")[1]
            file_list = nonempty_files(f"{filedir}/tNGS/LC*/LC*merge_count.xls")
            for file1 in file_list:
                if not has_header(file1):
                    continue
                dic_t={}
                sample_name=file1.split("/")[-1].split(".")[0]
                if sample_name not in dic_sp:
                    warn(f"skip sample without sp_ap_info entry: {sample_name}")
                    continue
                with open(file1, "r") as FILE:
                    try:
                        next(FILE)
                    except StopIteration:
                        warn(f"skip headerless file: {file1}")
                        continue
                    for line in FILE:
                        lines = line.strip("\n").split("\t")
                        a_id=lines[0]
                        dp=lines[1]
                        if a_id in list_id:
                            dic_t[a_id]=dp
                t_result=[pici,m_id,sample_name]
                list_ap1=dic_sp[sample_name].split("|")[0].split(",")
                list_ap2=dic_sp[sample_name].split("|")[1].split(",")
                for id in list_id:
                    if dic_t.get(id):
                        if id not in list_ap1:
                            fm+=float(dic_t[id])
                            if id not in list_ap2:
                                fz+=float(dic_t[id])
                        t_result.append(dic_t[id])
                    else:
                        t_result.append("0")
                text="\t".join(t_result)+"\n"
                RAW.write(text)
            radio_f=fz/fm if fm else 0
            if not fm:
                warn(f"{pici} has zero denominator for AP ratio")
            print(pici,m_id,radio_f)
    return


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-pos", required=True, help="need dir")
    parse.add_argument("-pos2", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.outfile, arg.outdir,arg.pos,arg.pos2)


if __name__ == "__main__":
    main()
