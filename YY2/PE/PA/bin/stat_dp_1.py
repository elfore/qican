import argparse
import glob
import json
import os


def extract_sample_name(file_path):
    stem = os.path.basename(file_path).split(".")[0]
    return "-".join(stem.split("_")[1].split("-")[1:])


def get_total_bases(json_file):
    with open(json_file, "r") as handle:
        data = json.load(handle)
    return data["summary"]["after_filtering"]["total_bases"]


def collect_files(patterns):
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches
    return []


def stat(infile, pos, outfile, outdir):
    dicinfo = {}
    with open(infile, "r") as handle:
        for line in handle:
            lines = line.strip("\n").split("\t")
            dicinfo[lines[2]] = lines[0] + "|" + lines[1]

    ordered_pos = []
    head = ["runID", "mechineID", "sampleID"]
    with open(pos, "r") as pos_handle:
        next(pos_handle)
        for line in pos_handle:
            lines = line.strip("\n").split("\t")
            ordered_pos.append(lines[0] + "|" + lines[1])
            head.append("-".join(lines[0:4]))

    with open(f"{outdir}/{outfile}", "w") as raw:
        raw.write("\t".join(head) + "\n")
        for filedir, run_info in dicinfo.items():
            pici = run_info.split("|")[0]
            total_bases = {}
            json_files = collect_files([
                f"{filedir}/PA/result/json/*.json",
                f"{filedir}/PA/result2/json/*.json",
                f"{filedir}/PA/work/*/*/*.json",
            ])
            depth_files = collect_files([
                f"{filedir}/PA/result/filter_base_depth/*base_depth.txt",
                f"{filedir}/PA/result2/filter_base_depth/*base_depth.txt",
                f"{filedir}/PA/work/*/*/*.trim_primer.base_depth.txt",
            ])
            for json_file in json_files:
                if json_file.endswith("fastp.json"):
                    continue
                sample_name = extract_sample_name(json_file)
                try:
                    total_bases[sample_name] = float(get_total_bases(json_file))
                except (KeyError, ValueError, json.JSONDecodeError) as exc:
                    print(f"WARNING: failed to parse total_bases from {json_file}: {exc}")
            for depth_file in depth_files:
                sample_name = extract_sample_name(depth_file)
                if sample_name not in total_bases:
                    print(f"WARNING: missing json total_bases for {sample_name}, skip {depth_file}")
                    continue
                pos_depth = {}
                with open(depth_file, "r") as handle:
                    for line in handle:
                        lines = line.strip("\n").split("\t")
                        pos_depth[lines[0] + "|" + lines[1]] = lines
                row = [pici, pici.split("_")[1], sample_name]
                scale = total_bases[sample_name]
                for pos_key in ordered_pos:
                    if pos_key in pos_depth:
                        row.append((float(pos_depth[pos_key][2]) / scale) * 10000)
                    else:
                        row.append(0)
                raw.write("\t".join(map(str, row)) + "\n")


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-pos", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.pos, arg.outfile, arg.outdir)


if __name__ == "__main__":
    main()
