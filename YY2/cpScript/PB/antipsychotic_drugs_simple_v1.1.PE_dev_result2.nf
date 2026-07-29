params.script_path = "/mnt/gpfs/Dataset/04.project/Bioinfo/22.antipsychotic_drugs/00.script"
params.database = " /mnt/gpfs/Dataset/04.project/Bioinfo/22.antipsychotic_drugs/02.database"
params.software = "/mnt/gpfs/Dataset/04.project/Bioinfo/22.antipsychotic_drugs/05.software"
params.version = "v1.1.0"

params.ref = "${params.database}/masked_fa/mask_homo.fa"
params.amplicon_bed = "${params.database}/primer_bed/v30/pharm_sort.v30.bed"
params.merge_bed = "${params.database}/primer_bed/v30/pharm_sort_merge.v30.bed"
params.primer_fa = "${params.database}/primer_bed/v30/pharm_v30.fa"

params.hotspot_bed = "${params.database}/final_depth_pos.bed.clinical"
params.rs2primer = "${params.database}/rs2seq_primer/drug_v2_rs2primer.tsv.clinical"
params.filter_mut = "${params.database}/homo_pos.tsv"
params.pos_af = "${params.database}/mut_thousheld/pos_af.txt.clinical"
params.pos_af_mapq0 = "${params.database}/mut_thousheld/pos_af_mapq0.txt.clinical"
params.rs_list = "${params.database}/mut_thousheld/rs.list"
params.mut_pos = "${params.database}/mut_thousheld/mut_pos.bed"
params.mut_pos_mapq0 = "${params.database}/mut_thousheld/mut_pos_mapq0.bed"


params.cyp_region = "${params.database}/cyp_genotying/gene_region.txt"
params.cyp_allele = "${params.database}/cyp_genotying/target_allele.txt"
params.cyp_info = "${params.database}/cyp_genotying/metabolic_type.txt"

params.hla_bed = "${params.database}/hla_final.bed"
params.hla_reform = "${params.database}/hla_type_reform_new.tsv"
params.hla_add_grep1 = "${params.database}/hla_grep/hla_add_r1.list"
params.hla_add_grep2 = "${params.database}/hla_grep/hla_add_r2.list"
params.hla_speed_grep1 = "${params.database}/hla_grep/need.seq1"
params.hla_speed_grep2 = "${params.database}/hla_grep/need.seq2"

dataset = []
params.sample.each{it ->  dataset.add(tuple(it.key, it.value.tokenize(",").get(0), it.value.tokenize(",").get(1)))}
datasets = Channel.from(dataset)


workflow{
    fastp_PE(datasets)

    bwa_PE(fastp_PE.out.map {sample, fq1, fq2, json ->  tuple(sample,fq1,fq2)})
}




process fastp_PE {
    publishDir "${params.outpath}/QC", mode: 'copy', pattern: "*.json"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}.clean.R1.fastq"), path("${sample_name}.clean.R2.fastq"), path("${sample_name}.json")
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/fastp --adapter_sequence TGGAATTCTCGGGTGCCAAGGAACTCCAGT --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAG -l 50 --cut_tail --cut_tail_window_size 1  --cut_tail_mean_quality 15  --cut_right --cut_right_window_size 5 --cut_right_mean_quality 20 -i ${fq1} -I ${fq2} -o ${sample_name}.clean.R1.fastq -O ${sample_name}.clean.R2.fastq -j ${sample_name}.json -h ${sample_name}.html
    """
}


process bwa_PE {
    publishDir "${params.outpath}/bam", mode: 'copy', pattern: "*.base_depth.txt"
    publishDir "${params.outpath}/bamfile", mode: 'copy', pattern: "*.sort.bam*"
    input:
    tuple val(sample_name), path("${sample_name}.clean.R1.fastq"), path("${sample_name}.clean.R2.fastq")
    output:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"), path("${sample_name}.base_depth.txt")
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:cygnus\\tSM:${sample_name}\\tPU:unit1' -t 6 ${params.ref} ${sample_name}.clean.R1.fastq ${sample_name}.clean.R2.fastq | /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools view -Sb - > ${sample_name}.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools sort ${sample_name}.bam -o ${sample_name}.sort.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools index ${sample_name}.sort.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/trim_primer_pair.py ${sample_name}.sort.bam ${params.primer_fa} ${params.amplicon_bed} ${sample_name}.trim_primer.bam ${params.filter_mut}
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools sort ${sample_name}.trim_primer.bam -o ${sample_name}.trim_primer.sort.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools index ${sample_name}.trim_primer.sort.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools depth -d 100000 -q 30 -aa -b ${params.merge_bed} ${sample_name}.sort.bam  > ${sample_name}.base_depth.txt
    """
}


















