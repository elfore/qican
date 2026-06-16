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
    seqtk_downsample(datasets)
    fastp_PE(seqtk_downsample.out)

    bwa_PE(fastp_PE.out.map {sample, fq1, fq2, json ->  tuple(sample,fq1,fq2)})
    QC_stat(fastp_PE.out.map {sample, fq1, fq2, json ->  tuple(sample,json)})
    
    QC_stat.out.map {sample, qc_stat ->  tuple(sample, qc_stat)}.set {QC_stat_out}
    bwa_PE.out.map {sample, bam, bai ->  tuple(sample, bam, bai)}.set {bam_PE_out}

    filter_bam(bam_PE_out)

    freebayes_call_trim(filter_bam.out, params.merge_bed, params.ref)
    // add_cut_mut(freebayes_call_trim.out)

    freebayes_call_trim_mapq0(filter_bam.out, params.merge_bed, params.ref)
    freebayes_call_trim.out.combine(freebayes_call_trim_mapq0.out, by:0).map {sample, map20_vcf, bam, bai, mapq0_vcf -> tuple(sample, map20_vcf, bam, bai, mapq0_vcf)}| set {cyp_analysis_in}
    
    count_50ins(bam_PE_out)
    cal_cyp2a6_type46(bam_PE_out)

    cyp_analysis_in.combine(cal_cyp2a6_type46.out, by:0).map {sample, map20_vcf, bam, bai, mapq0_vcf, cyp_46type_tsv -> tuple(sample, map20_vcf, bam, bai, mapq0_vcf, cyp_46type_tsv)}| set{cyp_analysis_in_cyp46}
    cyp_analysis_new(cyp_analysis_in_cyp46)


    cyp_analysis_in.combine(count_50ins.out, by:0).map {sample, map20_vcf, bam, bai, mapq0_vcf, count_50ins_tsv -> tuple(sample, map20_vcf, bam, bai, mapq0_vcf, count_50ins_tsv)}| set {mut_combine_50ins}
    mut_combine_50ins.combine(cal_cyp2a6_type46.out, by:0).map {sample, map20_vcf, bam, bai, mapq0_vcf, count_50ins_tsv, cyp_46type_tsv -> tuple(sample, map20_vcf, bam, bai, mapq0_vcf, count_50ins_tsv, cyp_46type_tsv)}| set {mut_combine_50ins_cyp46}
    combine_mut(mut_combine_50ins_cyp46)

    grep_hla_fq(seqtk_downsample.out)
    cutadapter_cut(grep_hla_fq.out)
    hla_analysis_new2(cutadapter_cut.out)

    bam_PE_out.combine(QC_stat_out, by:0).map {sample, bam, bai, qc_stat -> tuple(sample, bam, bai, qc_stat)}| set {sample_qc_stat_bwa_tmp}
    sample_qc_stat_bwa_tmp.combine(filter_bam.out, by:0).map {sample, bam, bai, qc_stat, trim_bam, trim_bai -> tuple(sample, bam, bai, qc_stat, trim_bam, trim_bai)}| set {sample_qc_stat_bwa}
    sample_qc_stat_bwa.combine(hla_analysis_new2.out, by:0).map {sample, bam, bai, qc_stat, trim_bam, trim_bai, hla_result -> tuple(sample, bam, bai, qc_stat, trim_bam, trim_bai, hla_result)} | set {sample_qc_stat_bwa_addhla}
    sample_qc_stat_bwa_addhla.combine(combine_mut.out, by:0).map {sample, bam, bai, qc_stat, trim_bam, trim_bai, hla_result, combine_vcf -> tuple(sample, bam, bai, qc_stat, trim_bam, trim_bai, hla_result, combine_vcf)} | set {sample_qc_stat_bwa_addhla_vcf}
    bam_stat(sample_qc_stat_bwa_addhla_vcf, params.amplicon_bed, params.merge_bed)
    
    // tuple val(sample_name), path("${sample_name}.bam_stat.xls"), path("${sample_name}_qc_merge.csv"), path("${sample_name}_check_info.tsv"), path("${sample_name}.regions.bed.gz"), path("${sample_name}.reads_depth.txt"), path("${sample_name}_combine_rs.vcf")
    bam_stat.out.combine(cyp_analysis_new.out, by:0).map {sample, bam_stat, qc_merge, check_info, region_depth, reads_depth, combine_vcf, cyp_result -> tuple(sample, combine_vcf, cyp_result)} | set {mk_report_tmp}
    // combine_mut.out.combine(cyp_analysis_new.out, by:0).map {sample, combine_vcf, cyp_result -> tuple(sample, combine_vcf, cyp_result)} | set {mk_report_tmp}
    mk_report_tmp.combine(hla_analysis_new2.out, by:0).map {sample, combine_vcf, cyp_result, hla_result -> tuple(sample, combine_vcf, hla_result, cyp_result)} | set {mk_report_in}
    mk_report(mk_report_in)

}


process count_50ins{
    publishDir "${params.outpath}/count_50ins", pattern: "*_count.tsv"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    output:
    tuple val(sample_name), path("${sample_name}_count.tsv")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/python  ${params.script_path}/deal_long_ins7.py      ${sample_name}.sort.bam  cyg    >    ${sample_name}_count.tsv 
    """
}

process cal_cyp2a6_type46{
    publishDir "${params.outpath}/cal_cyp2a6_type46", pattern: "*_cyp.tsv"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    output:
    tuple val(sample_name), path("${sample_name}_cyp.tsv")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/python ${params.script_path}/deal_cyp_46type2.py   ${sample_name}.sort.bam   ${sample_name}   >    ${sample_name}_cyp.tsv 
    """
}


process combine_mut{
    input:
    tuple val(sample_name), path("${sample_name}_trim_primer.freebayes.vcf"), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai"), path("${sample_name}_trim_primer.freebayes.mapq0.vcf"), path("${sample_name}_count.tsv"), path("${sample_name}_cyp.tsv")
    output:
    tuple val(sample_name), path("${sample_name}_combine.vcf")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/python  ${params.script_path}/reform_mut_new2.py ${sample_name}_trim_primer.freebayes.vcf  ${sample_name}_trim_primer.freebayes.mapq0.vcf  ${params.pos_af} ${params.pos_af_mapq0}  > ${sample_name}_combine.vcf
    tail -1 ${sample_name}_count.tsv >> ${sample_name}_combine.vcf
    tail -1 ${sample_name}_cyp.tsv >> ${sample_name}_combine.vcf
    """
}


process seqtk_downsample {
    publishDir "${params.outpath}/downsample_fq", pattern: "*.gz"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}.downsample.R1.fq.gz"), path("${sample_name}.downsample.R2.fq.gz")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/seqtk sample -s 1000 ${fq1} ${params.target_reads_num} | /usr/bin/pigz -p 4 > ${sample_name}.downsample.R1.fq.gz
    /mnt/gpfs/Users/fanlei/miniconda3/bin/seqtk sample -s 1000 ${fq2} ${params.target_reads_num} | /usr/bin/pigz -p 4 > ${sample_name}.downsample.R2.fq.gz
    """
}


process fastp_PE {
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
    publishDir "${params.outpath}/bam", pattern: "*.sort.bam*"
    input:
    tuple val(sample_name), path("${sample_name}.clean.R1.fastq"), path("${sample_name}.clean.R2.fastq")
    output:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:cygnus\\tSM:${sample_name}\\tPU:unit1' -t 6 ${params.ref} ${sample_name}.clean.R1.fastq ${sample_name}.clean.R2.fastq | /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools view -Sb - > ${sample_name}.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools sort ${sample_name}.bam -o ${sample_name}.sort.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools index ${sample_name}.sort.bam
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
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/QC_stat.py ${sample_name}.json ${sample_name}
    """
}


process filter_bam {
    publishDir "${params.outpath}/filter_bam", pattern: "*.trim_primer.sort.bam*"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai")
    output:
    tuple val(sample_name), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai")
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/trim_primer_pair.py ${sample_name}.sort.bam ${params.primer_fa} ${params.amplicon_bed} ${sample_name}.trim_primer.bam ${params.filter_mut}
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools sort ${sample_name}.trim_primer.bam -o ${sample_name}.trim_primer.sort.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools index ${sample_name}.trim_primer.sort.bam
    """
}


process bam_stat{
    publishDir "${params.outpath}/bam_stat", pattern: "*.bam_stat.xls"
    publishDir "${params.outpath}/hotspot_depth", mode: 'copy', pattern: "*.regions.bed.gz"
    publishDir "${params.outpath}/depth", mode: 'copy', pattern: "*.reads_depth.txt"
    publishDir "${params.outpath}/qc_merge", mode: 'copy', pattern: "*qc_merge.csv"
    publishDir "${params.outpath}/check_info", pattern: "*_check_info.tsv"
    publishDir "${params.outpath}/combine_mut", mode: 'copy', pattern: "*_combine_rs.vcf"
    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"), path("${sample_name}.QC_stat.xls"), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai"), path("${sample_name}_hla_result.tsv"),  path("${sample_name}_combine.vcf")
    path amplicon_bed
    path merge_bed
    output:
    tuple val(sample_name), path("${sample_name}.bam_stat.xls"), path("${sample_name}_qc_merge.csv"), path("${sample_name}_check_info.tsv"), path("${sample_name}.regions.bed.gz"), path("${sample_name}.reads_depth.txt"), path("${sample_name}_combine_rs.vcf")
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools depth -d 100000 -q 30 -aa -b ${merge_bed} ${sample_name}.sort.bam  > ${sample_name}.base_depth.txt
    MOSDEPTH_PRECISION=5 /mnt/gpfs/Users/fanlei/miniconda3/bin/mosdepth  --fast-mode  -b ${params.hotspot_bed}   ${sample_name}  ${sample_name}.trim_primer.sort.bam
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/bedtools coverage -f 0.45 -a ${amplicon_bed} -b ${sample_name}.sort.bam > ${sample_name}.reads_depth.txt
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/samtools flagstat ${sample_name}.sort.bam > ${sample_name}.flagstat
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python  ${params.script_path}/bam_stat_speed.py -bam ${sample_name}.sort.bam -base_depth ${sample_name}.base_depth.txt -amplicon_depth ${sample_name}.regions.bed.gz -read_depth_file ${sample_name}.reads_depth.txt -flagstat ${sample_name}.flagstat -sample ${sample_name}
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/merge_qc_new5.py  ${sample_name}.QC_stat.xls ${sample_name}.bam_stat.xls ${params.rs2primer} ${sample_name}_hla_result.tsv ${sample_name}_combine.vcf ${sample_name}_qc_merge.csv ${sample_name}_check_info.tsv
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/mk_hotspot_vcf2.py  ${sample_name}_combine.vcf ${params.rs_list} ${sample_name}.regions.bed.gz > ${sample_name}_combine_rs.vcf
    """
}


process freebayes_call_trim{
    publishDir "${params.outpath}/freebayes_trim", pattern: "*_trim_primer.freebayes.vcf"
    input:
    tuple val(sample_name), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai")
    path merge_bed
    path ref
    output:
    tuple val(sample_name), path("${sample_name}_trim_primer.freebayes.vcf"), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/freebayes -C 5 -F 0.005 -m 20 -q 20 -t ${params.mut_pos} -f ${params.ref} ${sample_name}.trim_primer.sort.bam > ${sample_name}_trim_primer.freebayes.vcf.tmp
    /mnt/gpfs/Users/fanlei/miniconda3/bin/bcftools norm  -f ${params.ref} -m-any ${sample_name}_trim_primer.freebayes.vcf.tmp -o ${sample_name}_trim_primer.freebayes.vcf
    """
}


process freebayes_call_trim_mapq0{
    publishDir "${params.outpath}/freebayes_trim_mapq0", pattern: "*_trim_primer.freebayes.mapq0.vcf"
    input:
    tuple val(sample_name), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai")
    path merge_bed
    path ref
    output:
    tuple val(sample_name), path("${sample_name}_trim_primer.freebayes.mapq0.vcf")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/freebayes -C 5 -F 0.01 -m 0  -t ${params.mut_pos_mapq0} -f ${params.ref} ${sample_name}.trim_primer.sort.bam > ${sample_name}_trim_primer.freebayes.vcf.tmp
    /mnt/gpfs/Users/fanlei/miniconda3/bin/bcftools norm  -f ${params.ref} -m-any ${sample_name}_trim_primer.freebayes.vcf.tmp -o ${sample_name}_trim_primer.freebayes.mapq0.vcf
    """
}

//可能需要删除
process add_cut_mut{
    publishDir "${params.outpath}/final_mut", pattern: "*_final.vcf"
    input:
    tuple val(sample_name), path("${sample_name}_trim_primer.freebayes.vcf"), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai")
    output:
    tuple val(sample_name), path("${sample_name}_final.vcf")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/python ${params.script_path}/grep_target_mut.py  ${sample_name}.trim_primer.sort.bam ${sample_name}_trim_primer.freebayes.vcf > ${sample_name}_final.vcf
    """
}


process cyp_analysis_new{
    errorStrategy 'ignore'
    publishDir "${params.outpath}/cyp_analysis", mode: 'copy', pattern: "*.cyp_allele.txt"
    input:
    tuple val(sample_name), path("${sample_name}_trim_primer.freebayes.vcf"), path("${sample_name}.trim_primer.sort.bam"), path("${sample_name}.trim_primer.sort.bam.bai"), path("${sample_name}_trim_primer.freebayes.mapq0.vcf"), path("${sample_name}_cyp.tsv")
    output:
    tuple val(sample_name), path("${sample_name}.cyp_allele.txt")
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/get_cyp_gt.py -map0 ${sample_name}_trim_primer.freebayes.mapq0.vcf -map20 ${sample_name}_trim_primer.freebayes.vcf -hot ${params.cyp_allele} -sample ${sample_name} -outdir ./
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/cyp_typing.py -info ${params.cyp_info} -mut ./${sample_name}.cyp_hot_spot_mut.txt -region ${params.cyp_region} -allele ${params.cyp_allele} -sample ${sample_name} -outdir ./ -cyp2a6_46 ${sample_name}_cyp.tsv
    """
}


process mk_report{
    publishDir "${params.outpath}/final_report", pattern: "*.报告.docx"
    input:
    tuple val(sample_name), path(combine_vcf), path(hla_result), path(cyp_result)
    output:
    tuple val(sample_name), path("${sample_name}.报告.docx")
    script:
    """
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/make_mut_result_v4.py -vcf ${combine_vcf} -hla ${hla_result} -cyp ${cyp_result} -outfile ${sample_name}.mut_info.txt
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/make_sample_info_v2.py -input1 ${sample_name}.mut_info.txt -outfile ${sample_name}.info.txt
    /mnt/gpfs/Users/huyangzhirong/02.sofware/miniconda3/bin/python ${params.script_path}/auto_report_AT_v2.py -inputfile ${sample_name}.info.txt -outputfile ${sample_name}.报告.docx
    """
}


process grep_hla_fq {
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}_grephla1.fq"), path("${sample_name}_grephla2.fq")
    script:
    """
    /mnt/gpfs/Users/fanlei/miniconda3/bin/python  ${params.script_path}/grep_target_fq_speed.py  ${params.hla_bed} ${params.hla_add_grep1} ${params.hla_add_grep2} ${params.hla_speed_grep1} ${params.hla_speed_grep2} /mnt/gpfs/Users/fanlei/miniconda3/bin/seqkit ${fq1} ${fq2} tmp_fq 3 2000 5 ${sample_name}_grephla1.fq ${sample_name}_grephla2.fq 
    """
}


process cutadapter_cut {
    publishDir "${params.outpath}/cut_fq", pattern: "*.gz"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}.cut.R1.fq.gz"), path("${sample_name}.cut.R2.fq.gz")
    script:
    """
    /mnt/gpfs/Users/zhangsijia/software/miniconda3/bin/cutadapt  -u  ${params.cut_num} -o ${sample_name}.cut.R1.fq.gz   ${fq1}
    /mnt/gpfs/Users/zhangsijia/software/miniconda3/bin/cutadapt  -u  ${params.cut_num} -o ${sample_name}.cut.R2.fq.gz   ${fq2}
    """
}


process hla_analysis_new2{
    errorStrategy 'ignore'
    publishDir "${params.outpath}/hla_analysis_new2", mode: 'copy', pattern: "*_hla_result.tsv"
    input:
    tuple val(sample_name), path(fq1), path(fq2)
    output:
    tuple val(sample_name), path("${sample_name}_hla_result.tsv")
    script:
    """
    export PATH=/mnt/gpfs/Users/fengbinxiao/00.App/miniconda3/envs/py310/bin/:\$PATH
    if [ \$(stat -c%s -L "${fq1}") -lt 1024 ] || [ \$(stat -c%s -L "${fq2}") -lt 1024 ]; then
    echo -e "Gene\tTypes\texon2_Depth20_cov\texon3_Depth20_cov\tDepth20_cov\nHLA-A_Allele1\tUnknown\tUnknown\tUnknown\tUnknown\nHLA-A_Allele2\tUnknown\tUnknown\tUnknown\tUnknown\nHLA-B_Allele1\tUnknown\tUnknown\tUnknown\tUnknown\nHLA-B_Allele2\tUnknown\tUnknown\tUnknown\tUnknown" > ${sample_name}_hla_result.tsv
    else
    /mnt/gpfs/Users/fengbinxiao/00.App/miniconda3/envs/py310/bin/python /mnt/gpfs/Users/fanlei/project/22.antipsychotic_drugs/00.script/speedup_hla/OptiTypePipeline.py -i ${fq1} ${fq2} -d -o ${sample_name}_hla -p ${sample_name} -v
    /mnt/gpfs/Users/fengbinxiao/00.App/miniconda3/envs/py310/bin/python ${params.script_path}/reform_optitype_new_qc_fix.py ${sample_name}_hla/${sample_name}_result.tsv ${params.hla_reform} result_exon_info.csv result_depth_info.csv > ${sample_name}_hla_result.tsv
    fi
    """
}
