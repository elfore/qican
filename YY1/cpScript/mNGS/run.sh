if [[ "${1:-}" == "__qsub__" ]]; then
  shift
  rule="${1:-unknown}"
  threads="${2:-1}"
  mem_gb="${3:-5}"
  jobscript="${4:-}"

  if [[ -z "$jobscript" ]]; then
    echo "ERROR: missing Snakemake jobscript argument" >&2
    exit 2
  fi

  if [[ -z "$threads" || "$threads" == "None" ]]; then
    threads=1
  fi

  if [[ -z "$mem_gb" || "$mem_gb" == "None" ]]; then
    mem_gb=5
  fi

  case "$rule" in
    bwa)
      threads="${BWA_THREADS:-2}"
      mem_gb="${BWA_MEM_GB:-6}"
      ;;
    bracken)
      threads="${BRACKEN_THREADS:-1}"
      mem_gb="${BRACKEN_MEM_GB:-5}"
      ;;
    kraken2)
      threads="${KRAKEN2_THREADS:-1}"
      mem_gb="${KRAKEN2_MEM_GB:-8}"
      ;;
    remain_reads_stat)
      threads="${REMAIN_THREADS:-1}"
      mem_gb="${REMAIN_MEM_GB:-2}"
      ;;
  esac

  exec qsub -V -cwd -q all.q -l hostname=!node03 -l vf="${mem_gb}G" -pe smp "${threads}" "$jobscript"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONDA_SH="/mnt/gpfs/Users/wangning/software/Miniconda/bin/activate"
CONDA_ENV="/mnt/gpfs/Users/luoshizhi/miniconda3/envs/metagenome"
SNAKEFILE="/mnt/gpfs1/Users/wangning/pipeline/metagenome/Snakefile"
CONFIGFILE="${CONFIGFILE:-sample.yaml}"
JOBS="${JOBS:-4}"
DEFAULT_MEM_GB="${DEFAULT_MEM_GB:-5}"
LATENCY_WAIT="${LATENCY_WAIT:-60}"

source "$CONDA_SH" "$CONDA_ENV"

exec snakemake -p \
  --configfile "$CONFIGFILE" \
  --default-resources "mem_gb=${DEFAULT_MEM_GB}" \
  --cluster "$SCRIPT_DIR/run.sh __qsub__ {rule} {threads} {resources.mem_gb}" \
  -j "$JOBS" \
  --latency-wait "$LATENCY_WAIT" \
  -s "$SNAKEFILE"
