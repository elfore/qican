import argparse
import glob
import os


def first_existing_file(patterns):
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[0]
    return None


def collect_existing_files(patterns):
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches
    return []


def get_pos_stat(infile, depth_files):
    pos_map = {}
    results = []
    with open(infile, "r") as pos_handle:
        for line in pos_handle:
            lines = line.strip("\n").split("\t")
            pos_map.setdefault(lines[0], {}).setdefault(lines[1], []).append(lines[3])

    for depth_file in depth_files:
        sample_name = os.path.basename(depth_file).split(".")[0]
        a = b = c = d = e = f = matched = missing = d30 = d50 = 0
        with open(depth_file, "r") as handle:
            for line in handle:
                lines = line.strip("\n").split("\t")
                if pos_map.get(lines[0]) and pos_map[lines[0]].get(lines[1]):
                    matched += 1
                    depth = int(lines[2])
                    if depth > 0:
                        a += 1
                    else:
                        missing += len(pos_map[lines[0]][lines[1]])
                    if depth >= 20:
                        b += 1
                    if depth >= 30:
                        d30 += 1
                    if depth >= 50:
                        d50 += 1
                    if depth >= 100:
                        c += 1
                    if depth >= 200:
                        d += 1
                    if depth >= 500:
                        e += 1
                    if depth >= 1000:
                        f += 1
        if matched == 0:
            status = "无有效位点"
            text = [sample_name, "0", str(missing), "0%(0)", "0%(0)", "0%(0)", "0%(0)", "0%(0)", "0%(0)", "0%(0)", "0%(0)", status]
        else:
            status = "合格" if d30 == matched else "不合格"
            text = [
                sample_name,
                str(matched),
                str(missing),
                f"{round((a / matched) * 100, 2)}%({a})",
                f"{round((b / matched) * 100, 2)}%({b})",
                f"{round((d30 / matched) * 100, 2)}%({d30})",
                f"{round((d50 / matched) * 100, 2)}%({d50})",
                f"{round((c / matched) * 100, 2)}%({c})",
                f"{round((d / matched) * 100, 2)}%({d})",
                f"{round((e / matched) * 100, 2)}%({e})",
                f"{round((f / matched) * 100, 2)}%({f})",
                status,
            ]
        results.append("\t".join(text) + "\n")
    return results


def stat(infile, pos, outfile, outdir):
    dicinfo = {}
    with open(infile, "r") as handle:
        for line in handle:
            lines = line.strip("\n").split("\t")
            dicinfo[lines[2]] = lines[0] + "|" + lines[1]

    with open(f"{outdir}/{outfile}", "w") as raw:
        for filedir, run_info in dicinfo.items():
            pici = run_info.split("|")[0]
            sample_stats = {}
            qc_file = first_existing_file([
                f"{filedir}/PA/result/QC/*final_QC_stat.xls",
                f"{filedir}/PA/result2/QC/*final_QC_stat.xls",
            ])
            bam_file = first_existing_file([
                f"{filedir}/PA/result/bam_stat/*final_bam_stat.xls",
                f"{filedir}/PA/result2/bam_stat/*final_bam_stat.xls",
            ])
            depth_files = collect_existing_files([
                f"{filedir}/PA/result/base_depth/*base_depth.txt",
                f"{filedir}/PA/result2/base_depth/*base_depth.txt",
                f"{filedir}/PA/work/*/*/*base_depth.txt",
            ])
            if not qc_file or not bam_file:
                print(f"WARNING: missing QC/bam_stat result under {filedir}/PA, skip")
                continue
            pos_dp = get_pos_stat(pos, depth_files)
            with open(qc_file, "r") as qc_handle, open(bam_file, "r") as bam_handle:
                next(qc_handle)
                next(bam_handle)
                for line in qc_handle:
                    lines = line.strip("\n").split("\t")
                    key = "-".join(lines[0].split("-")[1:]) if pici != "PL2402201" else lines[0].split("_")[1]
                    sample_stats.setdefault(key, []).extend(lines[1:])
                for line in bam_handle:
                    lines = line.strip("\n").split("\t")
                    key = "-".join(lines[0].split("-")[1:]) if pici != "PL2402201" else lines[0].split("_")[1]
                    sample_stats.setdefault(key, []).extend(lines[1:])
                for line in pos_dp:
                    lines = line.strip("\n").split("\t")
                    key = "-".join(lines[0].split("-")[1:]) if pici != "PL2402201" else lines[0].split("_")[1]
                    sample_stats.setdefault(key, []).extend(lines[1:])
            for key, values in sample_stats.items():
                machine_id = "-" if pici == "PL2402201" else pici.split("_")[1]
                raw.write(pici + "\t" + machine_id + "\t" + key + "\t" + "\t".join(map(str, values)) + "\n")


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
