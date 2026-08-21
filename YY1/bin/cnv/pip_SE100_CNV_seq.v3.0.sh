#baseline常染色体不区分性别
#与pip_SE150.sh的区别是，使用R脚本进行GC校正
fastp="/mnt/gpfs/Users/wangwenping/software/miniconda3/envs/WES/bin/fastp"
python="/mnt/gpfs/Users/wangwenping/software/Venv/python3/bin/python"
bwa="/mnt/gpfs/Software/Bioinfo/bwa/bwa"
adapter="ACCATCTCGGAGGTTGTTCCAGCGAAGAGT"
samtools="/mnt/gpfs/Software/Bioinfo/samtools/bin/samtools"
SampleGender="/mnt/gpfs/Users/huanghuichang/00.software/miniconda/envs/ngs_bits/bin/SampleGender"
bedtools="/mnt/gpfs/Users/luoshizhi/miniconda3/bin/bedtools"
Rscript="/mnt/gpfs/Software/Bioinfo/Miniconda3/bin/Rscript"
seqtk="/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/seqtk"
#database
bin1="/mnt/gpfs/Users/wangwenping/project/PGTA/bed/hg19_100kb_filter_mq.bed"
bin2="/mnt/gpfs/Users/wangwenping/project/PGTA/bed/1Mb_filter.bin"
ref="/mnt/gpfs/Users/wangning/database/human/hg19.fa"
mapq="/mnt/gpfs/Users/wangwenping/project/PGTA/mappingability/hg19_100kb_filter_mappability.txt"
mapq2="/mnt/gpfs/Users/wangwenping/project/PGTA/mappingability/1Mb_filter3_mappability.txt"
#baseline="/mnt/gpfs/Users/wangwenping/project/PGTA/baseline/CNV_seq/bin/baesline_100kb_stat_SE100.txt"
#baseline="/mnt/gpfs/Users/wangwenping/project/PGTA/analysis-new/SKII12330/rm_positive_region.baseline.xls"
#baseline2="/mnt/gpfs/Users/wangwenping/project/PGTA/baseline/CNV_seq/bin/1mb_cv_stat_SE100.txt" 
#baseline="/mnt/gpfs/Users/wangwenping/project/PGTA/baseline/CNV_seq/bin/baesline_100kb_stat_SE100_v2.txt"
#baseline2="/mnt/gpfs/Users/wangwenping/project/PGTA/baseline/CNV_seq/bin/1mb_cv_stat_SE100_v2.txt"
baseline="/mnt/gpfs/Users/huyangzhirong/01.project/13.PGTA/database/baseline/CNVseq/SE75/baseline_100kb_stat_SE75_2.5M_filter.txt"
baseline2="/mnt/gpfs/Users/huyangzhirong/01.project/13.PGTA/database/baseline/CNVseq/SE75/1mb_cv_stat_SE75_2.5M.txt"

samid=$1
R1=$2
reads=$3

source /mnt/gpfs/Users/wangwenping/software/Venv/python3/bin/activate
mkdir $samid
#####fastp
$fastp \
    -i $R1 -o $samid/${samid}.clean.R1.fq \
    -j $samid/${samid}.json \
    -w 10 \
    -a $adapter \
    -l 30 \
    -e 15 \
    -y \
    -h $samid/fastp.html

$seqtk sample $samid/${samid}.clean.R1.fq $reads > $samid/${samid}.clean.sub.fq

#######bwa and sort  官方说比对前不需要对reads进行过滤
`$bwa mem \
    -t 6 -M \
    -R '@RG\tID:$samid\tSM:$samid\tLB:$samid\tPL:PGS\tCN:cygnus' \
    $ref $samid/${samid}.clean.sub.fq \
    | $samtools view -bS - \
    > $samid/${samid}.bam
`

`
$samtools sort -T ./  -@ 10 -o $samid/${samid}.sort.bam  $samid/${samid}.bam
`
$samtools index $samid/${samid}.sort.bam 
rm $samid/${samid}.bam

$python /mnt/gpfs/Users/wangwenping/pythonProject/pipline/PGTA/new_pipeline/filter_map_len.py -bam $samid/${samid}.sort.bam -out $samid/${samid}.filter.bam -mapq 20
$samtools index $samid/${samid}.filter.bam
#########gender
$SampleGender -in $samid/${samid}.sort.bam  -method xy -min_male 0.18 -max_female 0.18 > $samid/${samid}.gender.txt 

#############rawbincount
$bedtools map -a $bin1 -b $samid/${samid}.filter.bam  -c 10 -o count > $samid/${samid}.100kb.readscount.txt 
$bedtools map -a $bin2 -b $samid/${samid}.filter.bam  -c 10 -o count > $samid/${samid}.1Mb.readscount.txt 

#########GC correction by R
$Rscript /mnt/gpfs/Users/huyangzhirong/01.project/13.PGTA/script/cbs.copy.R $samid/${samid}.100kb.readscount.txt $samid/${samid}.100kb.gc_correct.txt $samid
$Rscript /mnt/gpfs/Users/huyangzhirong/01.project/13.PGTA/script/cbs.copy.R $samid/${samid}.1Mb.readscount.txt $samid/${samid}.1Mb.gc_correct.txt $samid
###############normalization
$python /mnt/gpfs/Users/wangwenping/pythonProject/pipline/PGTA/new_pipeline/baseline_nor_v2.py -b $baseline -i $samid/${samid}.100kb.gc_correct.txt -g $samid/${samid}.gender.txt -m $mapq -o $samid/${samid}.cn_nor.txt
$python /mnt/gpfs/Users/wangwenping/pythonProject/pipline/PGTA/new_pipeline/baseline_nor_v2.py -b $baseline2 -i $samid/${samid}.1Mb.gc_correct.txt -g $samid/${samid}.gender.txt -m $mapq2 -o $samid/${samid}.1mb_cn_nor.txt

###############call cnv
$Rscript /mnt/gpfs/Users/huyangzhirong/01.project/13.PGTA/script/call_cbs_SE100.R $samid/${samid}.cn_nor.txt $samid/${samid}.cn.txt $samid
###############cnv scan
$python /mnt/gpfs/Users/huyangzhirong/01.project/13.PGTA/script/cnv_scan_mosaic.py -i $samid/${samid}.cn.txt -g $samid/${samid}.gender.txt -o $samid/${samid}.positive.txt -c /mnt/gpfs/Users/wangwenping/database/humandb/cytoBand.txt -r $samid/${samid}.result.txt

##################extractQC
$python /mnt/gpfs/Users/wangwenping/pythonProject/pipline/PGTA/new_pipeline/get_QC_v2.py -i $samid/${samid}.json -b $samid/${samid}.sort.bam -g $samid/${samid}.gender.txt -c $samid/${samid}.1mb_cn_nor.txt -c1 $samid/${samid}.cn_nor.txt -s $samtools -p $samid/${samid}.result.txt -o $samid/${samid}.totalqc.txt -n $samid -f  $samid/${samid}.filter.bam

##############pgs plot
python /mnt/gpfs/Users/wangwenping/pythonProject/pipline/PGTA/new_pipeline/cn_filter.py -i $samid/${samid}.cn.txt -o $samid/${samid}.cnplot.txt -p $samid/${samid}.positive.txt -g $samid/${samid}.gender.txt
/mnt/gpfs/Users/yangjinxurong/software/miniconda3/envs/stats/bin/Rscript /mnt/gpfs1/Users/yangjinxurong/projects/qican/YY1/bin/cnv/pgs_plot.R  $samid/${samid}.cnplot.txt $samid/${samid}.nor.png $samid