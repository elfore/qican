for i in $(cat sample);do echo "bash /mnt/gpfs1/Users/yangjinxurong/projects/qican/YY1/bin/cnv/pip_SE100_CNV_seq.v3.0.sh  ${i} /mnt/me4_nfs/02.analysis/260629161618_B173_SKII15393-R2JC-3-YY1-260629161619/Primary/FASTQ/L001/${i}_L001_R1.fq.gz  2500000" > run_${i}.sh;done
for i in $(cat sample);do  qsub -cwd -q all.q -pe smp 10 run_${i}.sh;done

