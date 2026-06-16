import argparse
import glob
import os


def collect_varscan_files(filedir):
    patterns = [
        f"{filedir}/PA/result/varscan/*.varscan.vcf",
        f"{filedir}/PA/result2/varscan/*.varscan.vcf",
        f"{filedir}/PA/work/*/*/*varscan.vcf",
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches
    return []


def stat(infile, outfile, outdir, pos):
    dicinfo = {}
    seen = {}
    pos_map = {}
    ordered_pos = []
    with open(infile, "r") as handle:
        for line in handle:
            lines = line.strip("\n").split("\t")
            dicinfo[lines[2]] = lines[0] + "|" + lines[1]

    head = ["runID", "mechineID", "sampleID"]
    with open(pos, "r") as pos_handle:
        next(pos_handle)
        for line in pos_handle:
            lines = line.strip("\n").split("\t")
            key = "|".join(lines[0:4])
            pos_map[key] = lines[3]
            ordered_pos.append(key)
            head.append("-".join(lines[0:4]))

    with open(f"{outdir}/{outfile}", "w") as raw:
        raw.write("\t".join(head) + "\n")
        for filedir, run_info in dicinfo.items():
            pici = run_info.split("|")[0]
            file_list = collect_varscan_files(filedir)
            if not file_list:
                print(f"WARNING: no varscan vcf found under {filedir}/PA")
                continue
            for file1 in file_list:
                sample_map = {}
                sample_field = os.path.basename(file1).split("_")[1].split(".")[0]
                if pici == "PL2402201":
                    sample_name = sample_field.split("-")[0]
                    sample_name2 = sample_field
                else:
                    parts = sample_field.split("-")
                    if len(parts) < 2:
                        print(f"WARNING: unexpected sample name format in {file1}")
                        continue
                    sample_name = parts[1]
                    sample_name2 = "-".join(parts[1:])
                sample_id = pici + "_" + sample_name2
                if sample_id in seen.setdefault(sample_name, set()) or os.path.getsize(file1) == 0:
                    continue
                seen[sample_name].add(sample_id)
                with open(file1, "r") as vcf_handle:
                    for line in vcf_handle:
                        if line.startswith("chr") and line.strip():
                            lines = line.strip("\n").split("\t")
                            freq = lines[9].split(":")[6]
                            if float(freq.split("%")[0]) > 3:
                                key = lines[0] + "|" + lines[1] + "|" + lines[3] + "|" + lines[4]
                                if key in pos_map and float(freq.split("%")[0]) >= 1:
                                    sample_map[key] = freq
                machine_id = "-" if pici == "PL2402201" else pici.split("_")[1]
                row = [pici, machine_id, sample_name2]
                row.extend(sample_map.get(var_key, "0%") for var_key in ordered_pos)
                raw.write("\t".join(row) + "\n")


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-infile", required=True, help="need dir")
    parse.add_argument("-pos", required=True, help="need dir")
    parse.add_argument("-outfile", required=True, help="need dir")
    parse.add_argument("--outdir", default=os.getcwd())
    arg = parse.parse_args()
    os.chdir(arg.outdir)
    stat(arg.infile, arg.outfile, arg.outdir, arg.pos)


if __name__ == "__main__":
    main()
