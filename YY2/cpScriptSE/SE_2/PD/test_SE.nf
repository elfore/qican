params.ref = "/mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/genome/hg19.mask.fa"
params.amplicon_bed = "/mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/bed/final_bed/final.bed"
params.merge_bed = "/mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/bed/final_bed/final.sort.merge.bed"
params.software = "/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin"

// 修改数据集读取逻辑：只保留一个 fq 路径
dataset = []
params.sample.each{it ->  dataset.add(tuple(it.key, it.value.tokenize(",").get(0)))}
datasets = Channel.from(dataset)

workflow{
    fastp_SE(datasets)
    // 修改 map 逻辑，只传递 sample 和单个 fq
    bwa_SE(fastp_SE.out[0].map {sample, fq, json ->  tuple(sample, fq)}, params.ref)
    QC_stat(fastp_SE.out[0].map {sample, fq, json ->  tuple(sample, json)})
    
    bwa_SE.out.combine(QC_stat.out, by:0).map {sample, bam, bai, qc_stat -> tuple(sample, bam, bai, qc_stat)} | set {bam_stat_in}
    
    bam_stat(bam_stat_in, params.amplicon_bed, params.merge_bed)
}

process fastp_SE {
    publishDir "${params.outpath}/dimer", pattern: "*.dimer_count.xls"
    input:
    tuple val(sample_name), path(fq1) // 仅输入 fq1
    output:
    tuple val(sample_name), path("${sample_name}.clean.sub.fastq"), path("${sample_name}.json")
    path "${sample_name}.dimer_count.xls"
    script:
    """
    # 去掉 -I 和 -O 参数，仅保留单端输入输出
    ${params.software}/fastp --adapter_sequence ${params.adapter} -l 50 -i ${fq1} -o ${sample_name}.clean.fastq -j ${sample_name}.json -h ${sample_name}.html
    ${params.software}/fastp --adapter_sequence ${params.adapter} -l 0 -i ${fq1} -o ${sample_name}.clean.dimer.fastq
    
    # 这里的 DimerStat 脚本如果内部逻辑是处理单端，则保持不变，此处仅根据单端输入修改参数
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/DimerStat.py -r ${sample_name}.clean.dimer.fastq -p /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/work/SKII12327/fasta/v13.fasta -o ./ -pf ${sample_name}
    
    # 抽样 50w reads
    ${params.software}/seqtk sample ${sample_name}.clean.fastq 500000 > ${sample_name}.clean.sub.fastq
    """
}

process bwa_SE {
    publishDir "${params.outpath}/bam", pattern: "*.sort.bam*"
    input:
    tuple val(sample_name), path("${sample_name}.clean.sub.fastq")
    path ref
    output:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    script:
    """
    find /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/genome -name "hg19.mask.fa.*" | awk '{print "ln -s "\$0 " ."}'|sh
    # BWA 仅接受一个 FASTQ 文件作为参数
    ${params.software}/bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:cygnus\\tSM:${sample_name}\\tPU:unit1' -t 8 ${ref} ${sample_name}.clean.sub.fastq | ${params.software}/samtools view -Sb - > ${sample_name}.bam
    ${params.software}/samtools sort ${sample_name}.bam -o ${sample_name}.sort.bam
    ${params.software}/samtools index ${sample_name}.sort.bam
    """
}

process QC_stat{
    publishDir "${params.outpath}/QC"
    input:
    tuple val(sample_name), path("${sample_name}.json")
    output:
    tuple val(sample_name), path("${sample_name}.QC_stat.xls")
    script:
    """
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/QC_stat.py ${sample_name}.json ${sample_name}
    """
}

process bam_stat{
    publishDir "${params.outpath}/depth", pattern: "*.reads_depth.txt"
    publishDir "${params.outpath}/depth_MAPQ", pattern: "*.reads_depth_MAPQ.txt"
    publishDir "${params.outpath}/base_depth", pattern: "*.base_depth.txt"
    publishDir "${params.outpath}/bam_stat", pattern: "*.bam_stat.xls"
    
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"), path("${sample_name}.QC_stat.xls")
    path amplicon_bed
    path merge_bed
    output:
    path "${sample_name}.bam_stat.xls"
    path "${sample_name}.reads_depth.txt"
    path "${sample_name}.reads_depth_MAPQ.txt"
    path "${sample_name}.base_depth.txt"
    
    script:
    """
    ${params.software}/samtools depth -q 13 -d 100000 -aa -b ${merge_bed} ${sample_name}.sort.bam  > ${sample_name}.base_depth.txt
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/amplicon_depth_stat.py -bed ${amplicon_bed} -bam ${sample_name}.sort.bam -outfile ${sample_name}.reads_depth.txt -workers 4
    ${params.software}/python /mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/PE150/SKII15275_SE_1/PD/amplicon_depth_stat_MAPQ.py  -bed ${amplicon_bed} -bam ${sample_name}.sort.bam -outfile ${sample_name}.reads_depth_MAPQ.txt -workers 4
    ${params.software}/samtools flagstat ${sample_name}.sort.bam > ${sample_name}.flagstat
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/normalize_depth.py -i ${sample_name}.QC_stat.xls -d ${sample_name}.reads_depth.txt -s ${sample_name}
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/bam_stat.py -bam ${sample_name}.sort.bam -base_depth ${sample_name}.base_depth.txt -amplicon_depth ${sample_name}.reads_depth.txt -bed ${amplicon_bed} -flagstat ${sample_name}.flagstat -hotspot /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/core_snp/core_pos.txt -sample ${sample_name}
    """
}