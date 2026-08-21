mkdir /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/05.Qican6/dingzhi_final/PE150/T7/PL260618-04/CN002342

cat /mnt/gpfs/Dataset/03.standard_seqdata/RawData/BGISEQ/PL260618-04/*CN002342*UDI02/*_R1.fastq.gz /mnt/gpfs/Dataset/03.standard_seqdata/RawData/BGISEQ/PL260618-04/*CN002342*UDI62*/*_R1.fastq.gz > /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/05.Qican6/dingzhi_final/PE150/T7/PL260618-04/CN002342/CP26000298-PL260618-04-CN002342-S01-D04-L02-UDI62_mixed_R1.fastq.gz

cat /mnt/gpfs/Dataset/03.standard_seqdata/RawData/BGISEQ/PL260618-04/*CN002342*UDI02/*_R2.fastq.gz /mnt/gpfs/Dataset/03.standard_seqdata/RawData/BGISEQ/PL260618-04/*CN002342*UDI62*/*_R2.fastq.gz > /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/05.Qican6/dingzhi_final/PE150/T7/PL260618-04/CN002342/CP26000298-PL260618-04-CN002342-S01-D04-L02-UDI62_mixed_R2.fastq.gz

bash run.sh

/mnt/gpfs/Users/yangjinxurong/software/miniconda3/envs/stats/bin/python /mnt/gpfs1/Users/yangjinxurong/projects/qican/QC6/bin/stat3.py -indir /mnt/gpfs/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/05.Qican6/dingzhi_final -infile /mnt/gpfs1/Users/yangjinxurong/projects/qican/QC6/bin/aa.txt