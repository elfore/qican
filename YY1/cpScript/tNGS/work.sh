for i in `cat sample`;do echo "sh /mnt/gpfs/Dataset/04.project/Bioinfo/20.tNGS/pipeline/V2/pip_ECC1_SE.sh $i /mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/SE75/GeneMind_0713/FASTQ/${i}_L000_R1.fq.gz " > $i.sh;done
for i in $(cat sample);do  qsub -cwd -q all.q -l h=!node03 -pe smp 8 ${i}.sh;done
# mkdir scripts
# mv LC*.sh LC*.sh.* scripts
