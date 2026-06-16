import argparse
import glob
import os


def get_ty(freq_text):
    return "阴性" if float(freq_text.split("%")[0]) <= 15 else "阳性"


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
    ref_calls = {}
    pos_map = {}
    ref_dir = ""
    seen = {}
    ordered_pos = []
    with open(infile, "r") as handle:
        for line in handle:
            lines = line.strip("\n").split("\t")
            if lines[0] == "PL2402201":
                ref_dir = lines[2]
            else:
                dicinfo[lines[2]] = lines[0] + "|" + lines[1]

    with open(pos, "r") as pos_handle:
        next(pos_handle)
        for line in pos_handle:
            lines = line.strip("\n").split("\t")
            key = "|".join(lines[0:4])
            pos_map[key] = lines[3]
            ordered_pos.append(key)

    with open(f"{outdir}/{outfile}", "w") as raw:
        raw.write("\t".join(["runID", "mechineID", "sampleID", "阳性符合率", "阴性符合率"]) + "\n")
        ref_summary = {}
        ref_files = collect_varscan_files(ref_dir) if ref_dir else []
        for file1 in ref_files:
            sample_name2 = os.path.basename(file1).split("_")[1].split(".")[0]
            yang = 0
            ying = 0
            with open(file1, "r") as vcf_handle:
                for line in vcf_handle:
                    if line.startswith("chr") and line.strip():
                        lines = line.strip("\n").split("\t")
                        key = lines[0] + "|" + lines[1] + "|" + lines[3] + "|" + lines[4]
                        if key in pos_map:
                            call = get_ty(lines[9].split(":")[6])
                            ref_calls.setdefault(sample_name2, {})[key] = call
                            if call == "阴性":
                                ying += 1
                            else:
                                yang += 1
            ref_calls.setdefault(sample_name2, {})
            for key in ordered_pos:
                if key not in ref_calls[sample_name2]:
                    ref_calls[sample_name2][key] = "阴性"
                    ying += 1
            ref_summary[sample_name2] = {"ying": ying, "yang": yang}

        for filedir, run_info in dicinfo.items():
            pici = run_info.split("|")[0]
            file_list = collect_varscan_files(filedir)
            if not file_list:
                print(f"WARNING: no varscan vcf found under {filedir}/PA")
                continue
            for file1 in file_list:
                sample_field = os.path.basename(file1).split("_")[1].split(".")[0]
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
                if sample_name2 not in ref_calls or sample_name2 not in ref_summary:
                    print(f"WARNING: missing PL2402201 reference calls for {sample_name2}, skip")
                    continue
                ying = 0
                yang = 0
                sample_calls = {}
                with open(file1, "r") as vcf_handle:
                    for line in vcf_handle:
                        if line.startswith("chr") and line.strip():
                            lines = line.strip("\n").split("\t")
                            key = lines[0] + "|" + lines[1] + "|" + lines[3] + "|" + lines[4]
                            if key in pos_map:
                                call = get_ty(lines[9].split(":")[6])
                                sample_calls[key] = call
                                if call == ref_calls[sample_name2][key] == "阳性":
                                    yang += 1
                                if call == ref_calls[sample_name2][key] == "阴性":
                                    ying += 1
                for key in ordered_pos:
                    if key not in sample_calls and ref_calls[sample_name2][key] == "阴性":
                        ying += 1
                yang_total = ref_summary[sample_name2]["yang"]
                ying_total = ref_summary[sample_name2]["ying"]
                yang_result = yang / yang_total if yang_total else 0
                ying_result = ying / ying_total if ying_total else 0
                machine_id = pici.split("_")[1]
                raw.write("\t".join(map(str, [pici, machine_id, sample_name2, yang_result, ying_result])) + "\n")


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
