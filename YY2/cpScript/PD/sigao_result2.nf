params.ref = "/mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/genome/hg19.mask.fa"
params.amplicon_bed = "/mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/bed/final_bed/final.bed"
params.merge_bed = "/mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/bed/final_bed/final.sort.merge.bed"
params.target_bed = "/mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/bed/final_bed/target_region.bed"
params.software = "/mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin"


dataset = []
params.sample.each{it ->  dataset.add(tuple(it.key, it.value.tokenize(",").get(0), it.value.tokenize(",").get(1)))}
datasets = Channel.from(dataset)


workflow{
    fastp_PE(datasets)
}


process fastp_PE {
    publishDir "${params.outpath}/QC", mode: 'copy', pattern: "*.json"
    publishDir "${params.outpath}/bam", mode: 'copy', pattern: "*.base_depth.txt"
    publishDir "${params.outpath}/bamfile", mode: 'copy', pattern: "*.sort.bam*"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}.clean.R1.fastq"), path("${sample_name}.clean.R2.fastq"), path("${sample_name}.json"), path("${sample_name}.base_depth.txt")

    script:
    """
    ${params.software}/fastp --adapter_sequence TGGAATTCTCGGGTGCCAAGGAACTCCAGT --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAG -l 50 -i ${fq1} -I ${fq2} -o ${sample_name}.clean.R1.fastq -O ${sample_name}.clean.R2.fastq -j ${sample_name}.json -h ${sample_name}.html
    find /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/genome -name "hg19.mask.fa.*" | awk '{print "ln -s "\$0 " ."}'|sh
    ${params.software}/bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:cygnus\\tSM:${sample_name}\\tPU:unit1' -t 8 ${params.ref} ${sample_name}.clean.R1.fastq ${sample_name}.clean.R2.fastq | ${params.software}/samtools view -Sb - > ${sample_name}.bam
    ${params.software}/samtools sort ${sample_name}.bam -o ${sample_name}.sort.bam
    ${params.software}/samtools index ${sample_name}.sort.bam
    ${params.software}/samtools depth -q 13 -d 100000 -aa -b ${params.merge_bed} ${sample_name}.sort.bam  > ${sample_name}.base_depth.txt

    """
}

















