mkdir -p results scripts

for i in $(cat sample); do
    printf '%s\n' "/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bwa mem -t 6 /mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/SE75/SKII12595/lambdaSD/ref/lambda_49SD.fa /mnt/me4_nfs/02.analysis/260818164738_B174_SKII15442-ZC-YY1-260818164739/Primary/FASTQ/L000/${i}_L000_R1.fq.gz | /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools sort -@ 10 -o ../results/${i}.sort.bam -" > "scripts/run_${i}.sh"
    printf '%s\n' "/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools index ../results/${i}.sort.bam" >> "scripts/run_${i}.sh"
done

cd scripts
for i in run_*.sh; do
    qsub -cwd -q all.q -l h='!node03' -pe smp 10 "$i"
done
