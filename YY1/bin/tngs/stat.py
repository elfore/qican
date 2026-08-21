import argparse
import os
import sys
import numpy as np


def stat(infile, infile2, outfile):
    list_id=[]
    dic={}
    with open(infile,"r") as IN1,open(infile2,"r") as IN2,open(outfile,"w") as OUT:
        for line in IN1:
            if line.startswith("amp"):continue
            id=line.strip("\n")
            list_id.append(id)
        for line in IN2:
            if line.startswith("amp"):
                OUT.write(line)
            else:
                lines=line.strip("\n").split("\t")
                dic[lines[0]]=line
        for key in list_id:
            text=dic[key]
            OUT.write(text)
    return

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-infile2", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile,arg.infile2, arg.outfile)


if __name__ == "__main__":
    main()