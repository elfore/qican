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
    bwa_PE(fastp_PE.out[0].map {sample, fq1, fq2, json ->  tuple(sample,fq1,fq2)}, params.ref)
    QC_stat(fastp_PE.out[0].map {sample, fq1, fq2, json ->  tuple(sample,json)})
    bwa_PE.out.combine(QC_stat.out, by:0).map {sample, bam, bai, qc_stat -> tuple(sample, bam, bai, qc_stat)}| set {sample_qc_stat_bwa}
    bam_stat(sample_qc_stat_bwa, params.amplicon_bed, params.merge_bed)
    bwa_PE.out.combine(bam_stat.out[3], by:0).map {sample, bam, bai, depth -> tuple(sample, bam, bai, depth)}|set {snp_call_in}
    snp_call(snp_call_in, params.target_bed, params.ref)
    fastp_PE.out[0].map {sample, cat_R1,cat_R2, json ->  tuple(sample, cat_R1,cat_R2)}.set {fastp_PE_out}
    mutscan_call(fastp_PE_out)
    // off_target_stat(bwa_PE.out, params.merge_bed)
    mutscan_call.out.combine(snp_call.out[0], by:0).map {sample, mutscan, snp -> tuple(sample, mutscan, snp)} |set {combine_mut_in}
    combine_mut(combine_mut_in)
    snp_call.out[1].combine(bam_stat.out[3], by:0).map {sample, vcf, depth -> tuple(sample, vcf, depth)} |set {genotyping_in}
    genotyping(genotyping_in)
    grep_hla_fq_new(bwa_PE.out)
    hla_analysis_new2(grep_hla_fq_new.out)
    combine_mut.out.combine(genotyping.out[0], by:0).combine(hla_analysis_new2.out, by:0) |set {report_in}
    report(report_in)
    summary(QC_stat.out.map {sample, qc_stat -> qc_stat}.collect(),bam_stat.out[0].collect(), bam_stat.out[2].collect(), bam_stat.out[3].map {sample, depth -> depth}.collect())
}


process fastp_PE {
    publishDir "${params.outpath}/dimer", pattern: "*.dimer_count.xls"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}.clean.R1.sub.fastq"), path("${sample_name}.clean.R2.sub.fastq"), path("${sample_name}.json")
    path "${sample_name}.dimer_count.xls"
    script:
    """
    ${params.software}/fastp --adapter_sequence TGGAATTCTCGGGTGCCAAGGAACTCCAGT --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAG -l 50 -i ${fq1} -I ${fq2} -o ${sample_name}.clean.R1.fastq -O ${sample_name}.clean.R2.fastq -j ${sample_name}.json -h ${sample_name}.html
    ${params.software}/fastp --adapter_sequence TGGAATTCTCGGGTGCCAAGGAACTCCAGT --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAG -l 0 -i ${fq1} -I ${fq2} -o ${sample_name}.clean.R1.dimer.fastq -O ${sample_name}.clean.R2.dimer.fastq
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/DimerStat.py -r ${sample_name}.clean.R1.dimer.fastq -p /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/work/SKII12327/fasta/v13.fasta -o ./ -pf ${sample_name}
    ${params.software}/seqtk sample ${sample_name}.clean.R1.fastq 500000 > ${sample_name}.clean.R1.sub.fastq
    ${params.software}/seqtk sample ${sample_name}.clean.R2.fastq 500000 > ${sample_name}.clean.R2.sub.fastq
    """
}


process bwa_PE {
    publishDir "${params.outpath}/bam", pattern: "*.sort.bam*"
    input:
    tuple val(sample_name), path("${sample_name}.clean.R1.sub.fastq"), path("${sample_name}.clean.R2.sub.fastq")
    path ref
    output:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    script:
    """
    find /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/genome -name "hg19.mask.fa.*" | awk '{print "ln -s "\$0 " ."}'|sh
    ${params.software}/bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:cygnus\\tSM:${sample_name}\\tPU:unit1' -t 8 ${ref} ${sample_name}.clean.R1.sub.fastq ${sample_name}.clean.R2.sub.fastq | ${params.software}/samtools view -Sb - > ${sample_name}.bam
    ${params.software}/samtools sort ${sample_name}.bam -o ${sample_name}.sort.bam
    ${params.software}/samtools index ${sample_name}.sort.bam
    """
}


process snp_call{
    publishDir "${params.outpath}/freebayes", pattern: "*.hotspots.xls"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"), path("${sample_name}.base_depth.txt")
    path target_bed
    path ref
    output:
    tuple val(sample_name), path("${sample_name}.hotspots.xls")
    tuple val(sample_name), path("${sample_name}.freebayes.vcf")
    script:
    """
    find /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/genome -name "hg19.mask.fa.*" | awk '{print "ln -s "\$0 " ."}'|sh
    /mnt/gpfs/Users/fanlei/miniconda3/bin/freebayes --haplotype-length -1 --min-coverage 20 -C 5 -F 0.005 -m 20 -q 20 -t ${target_bed} -f ${ref} ${sample_name}.sort.bam > ${sample_name}.freebayes.vcf.tmp
    /mnt/gpfs/Users/fanlei/miniconda3/bin/bcftools norm -f ${ref} -m-any ${sample_name}.freebayes.vcf.tmp -o ${sample_name}.freebayes.vcf
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/get_hotspot.py  -vcf ${sample_name}.freebayes.vcf -hot /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/hotspot/SNP.hotspot.txt -prefix ${sample_name} -depth ${sample_name}.base_depth.txt
    """
}

process mutscan_call{
    publishDir "${params.outpath}/mutscan"
    input:
    tuple val(sample_name),path("${sample_name}.clean.R1.sub.fastq"), path("${sample_name}.clean.R2.sub.fastq")
    output:
    tuple val(sample_name),path("${sample_name}.mutscan.txt")
    script:
    """
    /mnt/gpfs/Users/caiyilun/software/mutscan/mutscan -1 ${sample_name}.clean.R1.sub.fastq -2 ${sample_name}.clean.R2.sub.fastq -m /mnt/gpfs/Users/caiyilun/test/build_mutscan/sg_test/hot_gene.h > ${sample_name}.mutscan.tmp.txt
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/build_mutscan/mut_stat_v4.py -hot /mnt/gpfs/Users/caiyilun/test/build_mutscan/sg_test/hot_gene.h -mo ${sample_name}.mutscan.tmp.txt -outfile ${sample_name}.mutscan.txt -FHan /mnt/gpfs/Users/caiyilun/test/build_mutscan/sg_test/test.tsv
    """
}

process bam_stat{
    publishDir "${params.outpath}/bam_stat", mode: 'copy', pattern: "*.bam_stat.xls"
    publishDir "${params.outpath}/depth", mode: 'copy', pattern: "*.reads_depth.txt"
    publishDir "${params.outpath}/normal_depth", mode: 'copy', pattern: "*.normalized_raw_reads.xls"
    publishDir "${params.outpath}/base_depth", mode: 'copy', pattern: "*.base_depth.txt"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"), path("${sample_name}.QC_stat.xls")
    path amplicon_bed
    path merge_bed
    output:
    path "${sample_name}.bam_stat.xls"
    path "${sample_name}.reads_depth.txt"
    path "${sample_name}.normalized_raw_reads.xls"
    tuple val(sample_name), path("${sample_name}.base_depth.txt")
    script:
    """
    ${params.software}/samtools depth -q 13 -d 100000 -aa -b ${merge_bed} ${sample_name}.sort.bam  > ${sample_name}.base_depth.txt
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/amplicon_depth_stat.py -bed ${amplicon_bed} -bam ${sample_name}.sort.bam -outfile ${sample_name}.reads_depth.txt -workers 4
    ${params.software}/samtools flagstat ${sample_name}.sort.bam > ${sample_name}.flagstat
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/normalize_depth.py -i ${sample_name}.QC_stat.xls -d ${sample_name}.reads_depth.txt -s ${sample_name}
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/bam_stat.py -bam ${sample_name}.sort.bam -base_depth ${sample_name}.base_depth.txt -amplicon_depth ${sample_name}.reads_depth.txt -bed ${amplicon_bed} -flagstat ${sample_name}.flagstat -hotspot /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/core_snp/core_pos.txt -sample ${sample_name}
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

process combine_mut{
    publishDir "${params.outpath}/mutation", mode: 'copy'
    input:
    tuple val(sample_name), path("${sample_name}.mutscan.txt"), path("${sample_name}.hotspots.xls")
    output:
    tuple val(sample_name), path("${sample_name}.final_mut.txt")
    script:
    """
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/combine_mut.py -mutscan ${sample_name}.mutscan.txt -snp ${sample_name}.hotspots.xls -hotspot /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/hotspot/SNP.hotspot.txt -prefix ${sample_name}
    """
}

process genotyping{
    publishDir "${params.outpath}/genotyping", mode: 'copy'
    input:
    tuple val(sample_name), path("${sample_name}.freebayes.vcf"), path("${sample_name}.base_depth.txt")
    output:
    tuple val(sample_name), path("${sample_name}.allele.txt")
    path "${sample_name}.genotype_hotspot_mut.txt"
    script:
    """
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/get_gt.py -vcf ${sample_name}.freebayes.vcf -hot /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/genotying/target_allele.txt -depth ${sample_name}.base_depth.txt -prefix ${sample_name}
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/sg_genotying.py -mut ${sample_name}.genotype_hotspot_mut.txt -region /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/genotying/gene_region.txt -allele /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/genotying/target_allele.txt -sample ${sample_name} -outdir ./ -info /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/genotying/metabolic_type.txt
    """
}

process grep_hla_fq_new {
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    output:
    tuple val(sample_name), path("${sample_name}_HLA_cut.R1.fq.gz"), path("${sample_name}_HLA_cut.R2.fq.gz")
    script:
    """
    ${params.software}/samtools view -h ${sample_name}.sort.bam chr6:28477797-33448354 | ${params.software}/samtools sort -n /dev/stdin -@ 2 -o /dev/stdout | /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/sigao_config/software/samtools fastq /dev/stdin -1 ${sample_name}_HLA.R1.fq.gz -2 ${sample_name}_HLA.R2.fq.gz -s /dev/null -@ 2
    ${params.software}/python /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/sigao_config/script/longrange_primer.py /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/sigao_config/database/hla_longrange.fa /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/sigao_config/database/hla_longrange_grep.seq tmp_dir ${params.software}/seqkit ${sample_name}_HLA.R1.fq.gz ${sample_name}_HLA.R2.fq.gz 3 4 ${sample_name}_HLA_cut.R1.fq.gz ${sample_name}_HLA_cut.R2.fq.gz
    """
}


process hla_analysis_new2{
    publishDir "${params.outpath}/hla_analysis_new2", mode: 'copy', pattern: "*_hla_result.tsv"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}_hla_result.tsv")
    script:
    """
    export PATH=/mnt/gpfs/Users/fengbinxiao/00.App/miniconda3/envs/py310/bin/:\$PATH
    if [ \$(stat -c%s -L "${fq1}") -lt 1024 ] || [ \$(stat -c%s -L "${fq2}") -lt 1024 ]; then
    echo -e "Gene\tTypes\tDepth20_cov\nHLA-A_Allele1\tUnknown\tUnknown\nHLA-A_Allele2\tUnknown\tUnknown\nHLA-B_Allele1\tUnknown\tUnknown\nHLA-B_Allele2\tUnknown\tUnknown\nHLA-C_Allele1\tUnknown\tUnknown\nHLA-C_Allele2\tUnknown\tUnknown" > ${sample_name}_hla_result.tsv
    else
    /mnt/gpfs/Users/fengbinxiao/00.App/miniconda3/envs/py310/bin/python /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/speedup_hla/OptiTypePipeline.py -i ${fq1} ${fq2} -d -o ${sample_name}_hla -p ${sample_name} -v
    /mnt/gpfs/Users/fengbinxiao/00.App/miniconda3/envs/py310/bin/python /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/sigao_config/script/reform_optitype_longrange.py ${sample_name}_hla/${sample_name}_result.tsv /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/sigao_config/database/hla_type_reform.tsv result_exon_info.csv result_depth_info.csv > ${sample_name}_hla_result.tsv
    fi
    """
}



process off_target_stat{
    publishDir "${params.outpath}/off_target"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    path merge_bed
    output:
    path "${sample_name}.off_target.primer_cal.txt"
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools bamtobed -i ${sample_name}.sort.bam > ${sample_name}.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools sort -i ${sample_name}.bed > ${sample_name}.sorted.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools merge -i ${sample_name}.sorted.bed > ${sample_name}.sorted.merge.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools intersect -v -a ${sample_name}.sorted.merge.bed -b ${merge_bed} > ${sample_name}.off_target.bed
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools coverage -a ${sample_name}.off_target.bed -b ${sample_name}.sort.bam > ${sample_name}.off_target_depth.txt
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools intersect -a ${sample_name}.sort.bam -b ${sample_name}.off_target.bed > ${sample_name}.off_target.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools index ${sample_name}.off_target.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools fasta ${sample_name}.off_target.bam > ${sample_name}.off_target.fa
    /mnt/gpfs/Users/wangning/software/ncbi-blast-2.9.0+/bin/blastn -query ${sample_name}.off_target.fa -db /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/work/SKII12327/fasta/v13.fasta -out ${sample_name}.m6.out -evalue 1000 -word_size 15 -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore stitle"
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/10.CHD/pipeline/script/off_target_primer3.py ${sample_name}.m6.out ${sample_name}.off_target.primer_cal.txt
    """
}


process report{
    publishDir "${params.outpath}/report"
    input:
    tuple val(sample_name), path("${sample_name}.final_mut.txt"), path("${sample_name}.allele.txt"), path("${sample_name}_hla_result.tsv")
    output:
    path "${sample_name}.sg.docx"
    path "${sample_name}.xnxg.docx"
    script:
    """
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/sg_report_test/result_info_make/make_mut_result.py  -vcf ${sample_name}.final_mut.txt -cyp ${sample_name}.allele.txt -hla ${sample_name}_hla_result.tsv -outfile ${sample_name}.mut_result.sg.txt
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/sg_report_test/result_info_make/make_sample_info.py -input1 ${sample_name}.mut_result.sg.txt -outfile ${sample_name}.sample_info.sg.txt
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/sg_report_test/auto_report_sg.py -inputfile ${sample_name}.sample_info.sg.txt -outputfile ${sample_name}.sg.docx
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/xnxg_report_test/result_info_make/make_mut_result.py -vcf ${sample_name}.final_mut.txt -cyp ${sample_name}.allele.txt -hla ${sample_name}_hla_result.tsv -outfile ${sample_name}.mut_result.xnxg.txt
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/xnxg_report_test/result_info_make/make_sample_info.py -input1 ${sample_name}.mut_result.xnxg.txt -outfile ${sample_name}.sample_info.xnxg.txt
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/xnxg_report_test/auto_report_xnxg.py -inputfile ${sample_name}.sample_info.xnxg.txt -outputfile ${sample_name}.xnxg.docx
    """
}

process summary{
    input:
    path(x)
    path(y)
    path(z)
    path(m)
    script:
    """
    ${params.software}/python /mnt/gpfs/Users/caiyilun/test/auto_run_nf/test/need_file/hot_spot_stat.py -infile /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/database/core_snp/core_pos.txt -indir ${params.outpath}/base_depth --outdir ${params.outpath}/base_depth -bed ${params.amplicon_bed}
    ${params.software}/python /mnt/gpfs/Users/huyangzhirong/01.project/12.sigao/pipeline/script/sum_stat.py -indir ${params.outpath} --outdir ${params.outpath}
    """
}


