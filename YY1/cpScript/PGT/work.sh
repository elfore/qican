for i in $(cat sample);do echo "bash /mnt/gpfs1/Users/yangjinxurong/projects/qican/YY1/bin/pgta/PGTA.SE100_reanalysis_base.v0_6_0.sh  ${i} /mnt/me4_nfs/02.analysis/260702171042_B170_SKII15404-R2JC-ZY-YY-260702171043/Primary/FASTQ/L000/${i}_L000_R1.fq.gz SKII15404" > run_${i}.sh;done
for i in $(cat sample);do  qsub -cwd -q all.q -pe smp 10 -l h='!node03' run_${i}.sh;done
