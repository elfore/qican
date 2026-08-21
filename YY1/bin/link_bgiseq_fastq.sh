#!/usr/bin/env bash
set -euo pipefail

src_root="/mnt/gpfs/Dataset/03.standard_seqdata/RawData/BGISEQ/YY1-400M-V1.0.0"
dst_dir="/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/SE75/BGISEQ_0707/FASTQ"

dir_prefix="Sample_JZ26152227-YY1-400M-V1.0.0-"
file_prefix="JZ26152227-YY1-400M-V1.0.0-"

mkdir -p "$dst_dir"

linked=0
skipped=0
failed=0

for sample_dir in "$src_root"/"${dir_prefix}"*; do
  [[ -d "$sample_dir" ]] || continue

  dirname="$(basename "$sample_dir")"
  sample="${dirname#"$dir_prefix"}"

  for read in R1 R2; do
    src="${sample_dir}/${file_prefix}${sample}_combined_${read}.fastq.gz"
    dst="${dst_dir}/${sample}_L000_${read}.fq.gz"

    if [[ ! -f "$src" ]]; then
      echo "WARN: source missing, skip: $src" >&2
      failed=$((failed + 1))
      continue
    fi

    if [[ -e "$dst" || -L "$dst" ]]; then
      if [[ -L "$dst" && "$(readlink "$dst")" == "$src" ]]; then
        echo "OK: already linked: $dst"
        skipped=$((skipped + 1))
        continue
      fi

      echo "ERROR: destination exists and differs, skip: $dst" >&2
      failed=$((failed + 1))
      continue
    fi

    ln -s "$src" "$dst"
    echo "OK: linked: $dst -> $src"
    linked=$((linked + 1))
  done
done

echo "Summary: linked=$linked skipped=$skipped failed=$failed"

if [[ "$failed" -gt 0 ]]; then
  exit 1
fi