params.fastp="/mnt/gpfs/Users/wangning/software/fastp/fastp"
params.adapter="TGGAATTCTCGGGTGCCAAGGAACT"
params.python="/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python"
params.dimer_fa="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/primer/V15/V15_pool.primers.fa"
params.perl="/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/perl"
params.bwa="/mnt/gpfs/Software/Bioinfo/bwa/bwa"
params.panel_region="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/genome/V15/H37Rv_mask_IS6110_IS1081.addNC.add_NTM.add_cospecies.add_NTM_genotype.fasta"
params.samtools="/mnt/gpfs/Software/Bioinfo/samtools/bin/samtools"
params.seqkit="/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/seqkit"
params.NTM_fa="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/genome/NTM/V2/NTM_V2.fa"
params.species_name="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline/database/NTM/NTM.list.v2.txt"
params.allbed="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/bed/V15/V15.bed"
params.NCbed="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/bed/V3/PC.amplicon.bed"
params.bed="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/bed/V15/V15.sort.merge.bed"
params.drug_resistance_mutation="/mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/database/mutation/V14/hotspot.txt"
params.bedtools="/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools"
params.varscan="/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/varscan"


dataset = []
params.sample.each{it ->  dataset.add(tuple(it.key, it.value))}
datasets = Channel.from(dataset)


workflow{
    fastp(datasets)
    fastp.out.map {sample, fq, dimer,adapter ->  tuple(sample, fq, dimer,adapter)}|set {fastp_out}
    mk_bam(fastp.out.map {sample, fq, dimer,adapter ->  tuple(sample,fq)})
    mk_bam.out[0].map {sample, raw_bam, sort_bam,sort_bam_bai,  flagstat -> tuple(sample, sort_bam,sort_bam_bai)} | set {bam_out1}
    filter_bam(bam_out1)
    filter_bam.out.map {sample, sort_filter_bam -> tuple(sample, sort_filter_bam)} | set {filter_bam_out}
    ntm_indentify(fastp.out.map {sample, fq, dimer,adapter ->  tuple(sample,fq)})
   // varscan(filter_bam_out)
    filter_bam_out.combine(bam_out1, by:0).map {sample,sort_filter_bam,sort_bam,sort_bam_bai -> tuple(sample,  sort_filter_bam,sort_bam)} | set {qc_DNA1_in}
    qc_DNA1(qc_DNA1_in)
     qc_DNA1.out[0].combine(filter_bam_out, by:0).map {sample, amplicon_QC, amplicon_depth, NC_depth, sort_filter_bam -> tuple(sample, amplicon_QC,amplicon_depth,sort_filter_bam)} | set {qc_DNA2_tmp_in}
    qc_DNA2_tmp_in.combine(fastp_out, by:0).combine(mk_bam.out[1], by:0).map {sample, amplicon_QC,amplicon_depth,sort_filter_bam,fq,dimer,adapter, mut_depth, depth_stat, bam_stat, idxstats -> tuple(sample, amplicon_QC,amplicon_depth,sort_filter_bam,adapter, depth_stat, bam_stat)} | set {qc_DNA2_in}
    qc_DNA2(qc_DNA2_in)
    summary(qc_DNA2.out.map {sample, DNA_QC, filter_amplicon_depth,amplicon_depth_all, normalized_raw_reads ->normalized_raw_reads }.collect(), mk_bam.out[1].map {sample, mut_depth, depth_stat, bam_stat, idxstats -> idxstats}.collect())
}

process fastp {
    publishDir "${params.outpath}/dimer_stat" , pattern: "*dimer_count.xls"
    publishDir "${params.outpath}/adapter_stat" , pattern: "*adapter.xls"
    input:
    tuple val(sample_name), path(fq)
    output:
    tuple val(sample_name), path("${sample_name}.clean.fq.gz"),path("${sample_name}.dimer_count.xls"),path("${sample_name}.adapter.xls")
    script:
    """
    ${params.fastp}  -i ${fq} -o  ${sample_name}.clean.fq.gz  -a ${params.adapter} -l 50 --cut_right --cut_right_mean_quality 10 --failed_out ${sample_name}.fail.fq.gz -j ${sample_name}.json -h ${sample_name}.html
    ${params.python} /mnt/gpfs/Dataset/04.project/Bioinfo/09.XKMTB/01.pipeline/TB_V3/script/DimerStat.py -r ${sample_name}.fail.fq.gz -p ${params.dimer_fa} -b 10 -o ./ -pf ${sample_name}
    ${params.perl} /mnt/gpfs/Dataset/04.project/Bioinfo/09.XKMTB/01.pipeline/TB_V3/script/stat_adapter.pl ${sample_name}.json ${sample_name}
    """
}

process mk_bam {
    publishDir "${params.outpath}/bam", pattern: "*sort.bam*"
    publishDir "${params.outpath}/mutation", pattern: "*txt"
    input:
    tuple val(sample_name), path("${sample_name}.clean.fq.gz")
    output:
    tuple val(sample_name), path("${sample_name}.raw.bam"), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"),  path("${sample_name}.sort.flagstat")
    tuple val(sample_name), path("${sample_name}.mutations_with_depth.txt"), path("${sample_name}.depth_stats.txt"), path("${sample_name}.bam_reads_stat.txt"), path("${sample_name}.idxstats_depth.txt")
    script:
    """
    echo "run"
    ${params.bwa} mem -t 4 -M -R  '@RG\\tID:lib_lane\\tPL:cygnus\\tLB:lib\\tSM:${sample_name}' ${params.panel_region} ${sample_name}.clean.fq.gz | ${params.samtools} view -@ 4 -Sb - > ${sample_name}.raw.bam
    ${params.samtools} sort -@ 4 ${sample_name}.raw.bam --output-fmt BAM -o ${sample_name}.sort.bam
    ${params.samtools} index ${sample_name}.sort.bam
    ${params.samtools} flagstat ${sample_name}.sort.bam > ${sample_name}.sort.flagstat
    python /mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/script/bam_reads_stat.2.py ${sample_name}.sort.bam -o ${sample_name}.bam_reads_stat.txt
    ${params.python} /mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/script/mut_cov.py  ${params.drug_resistance_mutation} ${sample_name}.sort.bam ${sample_name}
    ${params.samtools} idxstats ${sample_name}.sort.bam > ${sample_name}.idxstats_depth.txt
    """
}


process filter_bam {
    publishDir "${params.outpath}/filter_bam"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    output:
    tuple val(sample_name), path("${sample_name}.sort.filter.bam")
    script:
    """
 #   ${params.python} /mnt/gpfs/Users/caiyilun/test/TB_test/filt_bam/filter_bam_p2.py  -bed ${params.allbed} -bam ${sample_name}.sort.bam -outfile ${sample_name}.sort.filter.bam -outfile2 ${sample_name}.sort.fail.bam
 #   ${params.samtools} sort -o ${sample_name}.sort.filter.tmp.bam ${sample_name}.sort.filter.bam
 #   ${params.samtools} sort -o ${sample_name}.sort.fail.tmp.bam ${sample_name}.sort.fail.bam
    mv ${sample_name}.sort.bam ${sample_name}.sort.filter.bam
    mv ${sample_name}.sort.bam.bai ${sample_name}.sort.filter.bam.bai
 #   ${params.samtools} index ${sample_name}.sort.filter.bam
 #   ${params.samtools} index ${sample_name}.sort.fail.bam
    """
}

process ntm_indentify {
    publishDir "${params.outpath}/ntm"
    input:
    tuple val(sample_name), path("${sample_name}.clean.fq.gz")
    output:
    //tuple val(sample_name), path("${sample_name}.ntm.stat.txt"), path("${sample_name}.ntm.result.txt")
    tuple val(sample_name), path("${sample_name}.dp.stat.txt")
    script:
    """
    ${params.bwa} mem -t 4 -M -R  '@RG\\tID:lib_lane\\tPL:cygnus\\tLB:lib\\tSM:${sample_name}' ${params.NTM_fa} ${sample_name}.clean.fq.gz | ${params.samtools} view -@ 4 -Sb - > ${sample_name}.NTM.bam
    ${params.samtools} sort -@ 4 ${sample_name}.NTM.bam --output-fmt BAM -o ${sample_name}.NTM.sort.bam
    ${params.samtools} index ${sample_name}.NTM.sort.bam
    ${params.python} /mnt/gpfs/Users/huyangzhirong/01.project/01.TB/new_pipeline/script/ntm_filter.bam.py -i ${sample_name}.NTM.sort.bam -o ${sample_name}.NTM.sort.filter.bam
    ${params.samtools} idxstats  ${sample_name}.NTM.sort.filter.bam > ${sample_name}.dp.stat.txt
    #${params.python} /mnt/gpfs/Users/caiyilun/test/TB_test/ntm_stat/ntm_stat.py   -infile ${params.species_name}  -infile2 ${sample_name}.dp.stat.txt -outfile ${sample_name}.ntm.stat.txt -outfile2 ${sample_name}.ntm.result.txt
    """
}

process varscan {
    publishDir "${params.outpath}/varscan"
    input:
    tuple val(sample_name), path("${sample_name}.sort.filter.bam")
    output:
    tuple val(sample_name), path("${sample_name}.varscan.mutation.xls")
    script:
    """
    ${params.samtools} mpileup -d 20000 -l ${params.bed} -f ${params.panel_region} ${sample_name}.sort.filter.bam > ${sample_name}.filter.mpileup
    ${params.varscan} mpileup2cns ${sample_name}.filter.mpileup --min-reads2 4 --min-coverage 20 --output-vcf 1 --min-var-freq 0.01 --strand-filter 0 --variants > ${sample_name}.varscan.raw.vcf
    ${params.python} /mnt/gpfs/Dataset/04.project/Bioinfo/09.XKMTB/01.pipeline/TB_V3/script/get_mutation.py ${params.drug_resistance_mutation} ${sample_name}.varscan.raw.vcf ${sample_name}.varscan.mutation.xls
    """
}

process qc_DNA1 {
    publishDir "${params.outpath}/QC_stat"
    input:
    tuple val(sample_name), path("${sample_name}.sort.filter.bam"),path("${sample_name}.sort.bam")
    output:
    tuple val(sample_name),path("${sample_name}.amplicon.QC.xls"), path("${sample_name}.amplicon.depth.xls"), path("${sample_name}.NC.depth.xls")
    path("${sample_name}.off_target_depth.txt")
    script:
    """
    ${params.python} /mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/script/QC_amplicon_reads_g.py  -i ${sample_name}.sort.bam -r ${params.bed} -o ./ -s ${params.samtools} -p ${sample_name}
    #${params.python} /mnt/gpfs/Dataset/04.project/Bioinfo/09.XKMTB/01.pipeline/TB_V3/script/QC_coverage.py   -s ${sample_name}  -p ./ -l ${params.drug_resistance_mutation}
    ${params.bedtools} coverage -f 0.4 -a ${params.allbed} -b ${sample_name}.sort.bam |cut -f 1-5 > ${sample_name}.amplicon.depth.xls
    ${params.bedtools} coverage -f 0.4 -a ${params.NCbed} -b ${sample_name}.sort.bam |cut -f 1-5 > ${sample_name}.NC.depth.xls
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools bamtobed -i ${sample_name}.sort.bam > ${sample_name}.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools sort -i ${sample_name}.bed > ${sample_name}.sorted.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools merge -i ${sample_name}.sorted.bed > ${sample_name}.sorted.merge.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools intersect -v -a ${sample_name}.sorted.merge.bed -b ${params.bed} > ${sample_name}.off_target.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools coverage -a ${sample_name}.off_target.bed -b ${sample_name}.sort.bam > ${sample_name}.off_target_depth.txt
    """
}

process qc_DNA2 {
    publishDir "${params.outpath}/QC_stat"
    input:
    tuple val(sample_name), path("${sample_name}.amplicon.QC.xls"), path("${sample_name}.amplicon.depth.xls"),path("${sample_name}.sort.filter.bam"),path("${sample_name}.adapter.xls"),path("${sample_name}.depth_stats.txt"), path("${sample_name}.bam_reads_stat.txt")
    output:
    tuple val(sample_name), path("${sample_name}.DNA.Q20.QC.xls"), path("${sample_name}.filter.amplicon.depth.xls"), path("${sample_name}.amplicon.depth.all.xls"), path("${sample_name}.normalized_raw_reads.xls")
    script:
    """
    ${params.python}  /mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/script/DNA_BAM_QC.py  ${sample_name}.sort.filter.bam  > ${sample_name}.DNA.Q20.QC.xls
    ${params.bedtools} coverage -f 0.4  -a ${params.allbed} -b ${sample_name}.sort.filter.bam | cut -f 1-5> ${sample_name}.filter.amplicon.depth.xls
    ${params.python} /mnt/gpfs/Users/caiyilun/test/TB_test/amplicon_depth.py   -r ${sample_name}.amplicon.depth.xls -f ${sample_name}.filter.amplicon.depth.xls -o ./ -s ${sample_name}
    ${params.python} /mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/script/normalize.py  -i ${sample_name}.DNA.Q20.QC.xls -d ${sample_name}.amplicon.depth.all.xls -o ./ -s ${sample_name}
    """
}

process summary{
    input:
    path(x)
    path(y)
    script:
    """
    ${params.python}  /mnt/gpfs/Users/huyangzhirong/01.project/01.TB/pipeline_TB2/script/stat_g.py  -indir ${params.outpath} --outdir ${params.outpath}
    """
}
