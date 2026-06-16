params.ref = "/mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/genome/hg19/hg19.fa"
params.amplicon_bed = "/mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/bed/02.PE/V26/v26.amplicon.sort.bed"
params.merge_bed = "/mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/bed/02.PE/V26/v26.amplicon.sort.merge.bed"
params.software = "/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin"


dataset = []
params.sample.each{it ->  dataset.add(tuple(it.key, it.value.tokenize(",").get(0), it.value.tokenize(",").get(1)))}
datasets = Channel.from(dataset)


workflow{
    fastp_PE(datasets)
    bwa_PE(fastp_PE.out[0].map {sample, fq1, fq2, json ->  tuple(sample,fq1,fq2)}, params.ref)
    QC_stat(fastp_PE.out[0].map {sample, fq1, fq2, json ->  tuple(sample,json)})

    QC_stat.out.map {sample, qc_stat ->  tuple(sample, qc_stat)}.set {QC_stat_out}
    bwa_PE.out.map {sample, bam, bai ->  tuple(sample, bam, bai)}.set {bam_PE_out}
    bam_PE_out.combine(QC_stat_out, by:0).map {sample, bam, bai, qc_stat -> tuple(sample, bam, bai, qc_stat)}| set {sample_qc_stat_bwa}
    bam_stat(sample_qc_stat_bwa, params.amplicon_bed, params.merge_bed)
    filter_bam(bwa_PE.out, params.ref, params.amplicon_bed, params.merge_bed)
}


process fastp_PE {
    publishDir "${params.outpath}/json", mode: 'copy', pattern: "*.json"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}.clean.R1.fastq"), path("${sample_name}.clean.R2.fastq"), path("${sample_name}.json")
    script:
    """
    ${params.software}/fastp --adapter_sequence TGGAATTCTCGGGTGCCAAGGAACTCCAGT --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAG -l 50 -i ${fq1} -I ${fq2} -o ${sample_name}.clean.R1.fastq -O ${sample_name}.clean.R2.fastq -j ${sample_name}.json -h ${sample_name}.html
    """
}


process bwa_PE {
    publishDir "${params.outpath}/bam", pattern: "*.sort.bam*"
    input:
    tuple val(sample_name), path("${sample_name}.clean.R1.fastq"), path("${sample_name}.clean.R2.fastq")
    path ref
    output:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    script:
    """
    find /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/genome/hg19 -name "hg19.fa.*" | awk '{print "ln -s "\$0 " ."}'|sh
    ${params.software}/bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:cygnus\\tSM:${sample_name}\\tPU:unit1' -t 8 ${ref} ${sample_name}.clean.R1.fastq ${sample_name}.clean.R2.fastq | ${params.software}/samtools view -Sb - > ${sample_name}.bam
    ${params.software}/samtools sort ${sample_name}.bam -o ${sample_name}.sort.bam
    ${params.software}/samtools index ${sample_name}.sort.bam
    """
}


process bam_stat{
    publishDir "${params.outpath}/bam_stat", mode: 'copy', pattern: "*.bam_stat.xls"
    publishDir "${params.outpath}/depth", mode: 'copy', pattern: "*.reads_depth.txt"
    publishDir "${params.outpath}/normal_depth", pattern: "*.normalized_raw_reads.xls"
    publishDir "${params.outpath}/base_depth", mode: 'copy', pattern: "*.base_depth.txt"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"), path("${sample_name}.QC_stat.xls")
    path amplicon_bed
    path merge_bed
    output:
    path "${sample_name}.bam_stat.xls"
    path "${sample_name}.reads_depth.txt"
    path "${sample_name}.normalized_raw_reads.xls"
    path "${sample_name}.base_depth.txt"
    script:
    """
    ${params.software}/samtools depth -q 13 -d 100000 -aa -b ${merge_bed} ${sample_name}.sort.bam  > ${sample_name}.base_depth.txt
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/chd_20/ap_dep_stat/ap_dep_stat_p.py -bed ${amplicon_bed} -bam ${sample_name}.sort.bam -outfile ${sample_name}.reads_depth.txt
    ${params.software}/samtools flagstat ${sample_name}.sort.bam > ${sample_name}.flagstat
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/normalize_depth.py -i ${sample_name}.QC_stat.xls -d ${sample_name}.reads_depth.txt -s ${sample_name}
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/bam_stat.py -bam ${sample_name}.sort.bam -base_depth ${sample_name}.base_depth.txt -amplicon_depth ${sample_name}.reads_depth.txt -bed ${amplicon_bed} -flagstat ${sample_name}.flagstat -sample ${sample_name}
    """
}



process QC_stat{
    publishDir "${params.outpath}/QC", mode: 'copy'
    input:
    tuple val(sample_name), path("${sample_name}.json")
    output:
    tuple val(sample_name), path("${sample_name}.QC_stat.xls")
    script:
    """
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/QC_stat.py ${sample_name}.json ${sample_name}
    """
}

process filter_bam {
    publishDir "${params.outpath}/filter_bam", pattern: "*.trim_primer.sort.bam*"
    publishDir "${params.outpath}/filter_base_depth", mode: 'copy', pattern: "*.trim_primer.base_depth.txt"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    path ref
    path amplicon_bed
    path merge_bed
    output:
    tuple val(sample_name), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai")
    path "${sample_name}.trim_primer.base_depth.txt"
    script:
    """
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/trim_primer_pair.py ${sample_name}.sort.bam /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/work/02.PE/SKII12081/total_primers.fasta ${amplicon_bed} ${sample_name}.trim_primer.bam /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/test/PE/2.0/homo_pos.tsv
    ${params.software}/samtools sort ${sample_name}.trim_primer.bam -o ${sample_name}.trim_primer.sort.bam
    ${params.software}/samtools index ${sample_name}.trim_primer.sort.bam
    ${params.software}/samtools depth -q 13 -d 100000 -aa -b ${merge_bed} ${sample_name}.trim_primer.sort.bam  > ${sample_name}.trim_primer.base_depth.txt
    """
}

