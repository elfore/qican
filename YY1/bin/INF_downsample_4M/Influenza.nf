#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import java.io.File
import groovy.json.JsonSlurper


dataset = []
if (params.reffa){
    println ("")
    datasets = Channel.from(dataset)
}else {
    params.sample.each { sampleitem ->
        fq_list = sampleitem.fqpath.tokenize(",")
        def local_fq
        if (fq_list.size() == 1) {
            local_fq = fq_list.get(0)
        } else if (fq_list.size() == 2) {
            for(item in fq_list){
                if (item.contains("Read1")){
                    local_fq = item
                }
            }
        }
        dataset.add(tuple(sampleitem.sampleid, local_fq))
    }
    datasets = Channel.from(dataset)

    params.Runid = params.sample[0].runid
    params.sample_name = params.sample[0].sampleid
    params.resultdir = "${params.outpath}/Result/${params.sample_name}"
    params.reportdir = "${params.outpath}/Result"
    params.NA_type = params.sample[0].NA_type
    params.Seq_type = params.sample[0].seq_type
    params.Control_type = params.sample[0].control_type
}
params.thread = 10
status_file = file('analysis_status.txt')

workflow {
    status_file.text =  "processing\n"
    qualityControl(datasets)
    bwa(qualityControl.out[0].map{sample, read, json -> tuple(sample, read)}, params.genome)
    spades(bwa.out[0], params.genome, params.database)

    //target DB
    qualityControl.out[0].combine(spades.out[0], by:0).map {sample, read, json, fasta -> tuple(sample, read, json, fasta)} |set{to_bwa}
    bwa1(to_bwa, params.database)

    qualityControl.out[0].combine(bwa1.out[0],by:0).combine(spades.out[0], by:0).map{sample,read,json, bam, bai, qc, depth, refDB -> tuple(sample, bam, bai, qc, depth, refDB,json)} | set{call_input}
    cal_variation(call_input)

    cal_variation.out[0].combine(spades.out[0], by:0).map{ sample, vcf, IGV,consensus,fasta -> tuple(sample, vcf, fasta) } | set{ to_anno_mut }
    anno_mut(to_anno_mut, params.genome)

    bwa1.out[0].combine(cal_variation.out[0], by:0).map{sample, bam, bai, qc, depth, vcf,IGV, consensus -> tuple(sample, qc, consensus)} | set{call_phylo}
    phylogenetic(call_phylo, params.database)

    //small DB
    small_out = small(to_bwa, params.database)
    report1_out = report1(
        small.out[0].map{sample, qc, chromQC, summaryQC, fastpQC, result, png, genotype-> chromQC}.collect(),
        small.out[0].map{sample, qc, chromQC, summaryQC, fastpQC, result, png, genotype-> summaryQC}.collect(),
        small.out[0].map{sample, qc, chromQC, summaryQC, fastpQC, result, png, genotype-> fastpQC}.collect(),
        small.out[0].map{sample, qc, chromQC, summaryQC, fastpQC, result, png, genotype-> result}.collect()
    )

    guaranteed_ch = report1_out.ifEmpty ( Channel.of('done') )

    report_out = report(
        bwa1.out[1].map{sample, chromQC,png -> chromQC }.collect(),
        cal_variation.out[1].map{sample, summaryQC, result,  variantstat, fastpQC -> summaryQC }.collect(),
        cal_variation.out[1].map{sample, summaryQC, result,  variantstat, fastpQC -> fastpQC }.collect(),
        cal_variation.out[1].map{sample, summaryQC, result,  variantstat, fastpQC -> result }.collect(),
        cal_variation.out[1].map{sample, summaryQC, result,  variantstat, fastpQC -> variantstat }.collect(),
        anno_mut.out[0].map{sample, anno -> anno}.collect(),
        guaranteed_ch
    )

    report_ch = report_out.ifEmpty ( Channel.of('done') )
    json_result(report_ch, guaranteed_ch)

    if (params.reffa){
        phylogenetic_fa(params.reffa, params.database)
    }
}


workflow.onError {
    status_file.text = "${ workflow.success ? 'done\n' : 'error\n' }"
}
workflow.onComplete = {
    status_file.text = "${ workflow.success ? 'done\n' : 'error\n' }"
}

// ==============================================
process qualityControl {
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.json", mode: 'copy'
//    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.clean.fq.gz", mode: 'symlink'

    input:
    tuple val(sample_name), path(read)

    output:
    tuple val(sample_name), path("${sample_name}.clean.fq.gz"), path("${sample_name}.json")

    script:
    """
    ${params.fastp} \
        -i $read \
        -o ${sample_name}.clean.fq.gz \
        -h ${sample_name}.fastp.html \
        -j ${sample_name}.json -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --reads_to_process 4000000 \
        -w $params.thread
    """
}

// ==============================================
process bwa {
    input:
    tuple val(sample_name), path("${sample_name}.clean.fq.gz")
    val params.genome

    output:
    tuple val(sample_name), path("${sample_name}.mapping.fq")

    script:
    """
    ${params.bwa} mem -t $params.thread  \
        ${params.genome}/largeDB/Influenza_wholeDB.fasta \
        ${sample_name}.clean.fq.gz > ${sample_name}.sam

    ${params.samtools} view  -@ $params.thread -bS ${sample_name}.sam > ${sample_name}.bam
    ${params.samtools} sort -@ $params.thread  ${sample_name}.bam > ${sample_name}.largeDB.bam
    ${params.samtools} index -@ $params.thread ${sample_name}.largeDB.bam
    ${params.samtools} bam2fq -F 0x4  ${sample_name}.largeDB.bam > ${sample_name}.mapping.fq
    """
}

// ==============================================
process spades {
//    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*fasta", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*xls", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.mapping.fq")
    val params.genome
    val params.database

    output:
    tuple val(sample_name), path("refDB.fasta")
    tuple val(sample_name), path("blast_ref_stat.xls"),path("blast_ref_stat.scored.xls")

    script:
    """
   ${params.spades_python} ${params.spades_py} -s ${sample_name}.mapping.fq \
     -o spades/ -t 5 --rnaviral

    if [ -s "spades/scaffolds.fasta" ];then
        ${params.blastn} \
            -query spades/scaffolds.fasta -db ${params.genome}/largeDB/Influenza_wholeDB.fasta \
            -out blast.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qlen qstart qend sstrand slen sstart send evalue bitscore stitle" \
            -evalue 0.1  -max_target_seqs 50000 -perc_identity 90  -dust no -num_threads $params.thread

        ${params.perl} ${params.script}/filter_blast.pl \
            -i blast.txt \
            -o ./ \
            -r ${params.genome}/largeDB/Influenza_wholeDB.fasta \
            -g ${params.database}/GCF_seg_info.txt \
            -s ${params.database}/segment.txt \
            -t ${sample_name}
    fi

    if [ ! -f "refDB.fasta" ];then
        touch refDB.fasta
    else
       ${params.python} ${params.script}/flu_assembly_segs.py \\
           --segdedup-complete  ${params.database}/GCF_seg_info.txt \\
           --dedup-complete  ${params.database}/dedup.GCF_seg_info.txt \\
           --assemblies assemble.id \\
           -o  map.tsv
      ${params.python} ${params.script}/flu_restore_id.py --map map.tsv -i refDB.fasta -o refDB.fasta.tmp
      mv refDB.fasta.tmp refDB.fasta

    fi
    if [ ! -f "blast_ref_stat.xls" ];then
        touch blast_ref_stat.xls
    fi
    if [ ! -f "blast_ref_stat.scored.xls" ];then
       touch blast_ref_stat.scored.xls
    fi
    """
}

// ==============================================
process small {
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.fastpQC.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.result.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.summaryQC.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.chromQC.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.qc.xls", mode: 'copy'

    input:
    tuple val(sample_name), path(clean_fq), path(json), path(refDB)
    val params.database

    output:
    tuple val(sample_name), path("${sample_name}.qc.xls"),path("${sample_name}.chromQC.xls"), path("${sample_name}.summaryQC.xls"), path("${sample_name}.fastpQC.xls"), path("${sample_name}.result.xls"), path("${sample_name}.cov.png"),path("${sample_name}.genotype.tsv")

    when:
     refDB.size() == 0

    script:
    """
    ${params.bwa} mem -t $params.thread -M \
        ${params.database}/blastDB/blastDB.fasta \
        ${sample_name}.clean.fq.gz > ${sample_name}.sam

    ${params.samtools} view -@ $params.thread -Sb ${sample_name}.sam | \
    ${params.samtools} sort -@ $params.thread - > ${sample_name}.sort.bam

    ${params.samtools} index ${sample_name}.sort.bam

    ${params.samtools} depth \
        -a -aa -b ${params.database}/blastDB/blastDB.bed \
        ${sample_name}.sort.bam > map.depth.txt

    ${params.python} ${params.script}/quality_control_v4.py \
        -bam ${sample_name}.sort.bam \
        -json ${sample_name}.json \
        -out ${sample_name}.qc.xls \
        -depth map.depth.txt \
        -ref ${params.database}/blastDB/blastDB.fasta \
        -chrom ${sample_name}.chromQC.xls

    ${params.Rscript_cov} ${params.script}/genome_coverage_DB.R \
        ${sample_name}.qc.xls map.depth.txt ${sample_name}.cov.png

    echo -e "seqName\tclade\tqc.overallStatus\tcoverage" > ${sample_name}.genotype.tsv

    ${params.python} ${params.script}/trans_result.py \
        -qc ${sample_name}.qc.xls \
        -json ${sample_name}.json \
        -genotype ${sample_name}.genotype.tsv \
        -outdir ./ \
        -sample $sample_name
    """
}

// ==============================================
process bwa1 {
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.html", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.cov.png", mode: 'copy'

    input:
    tuple val(sample_name), path(clean_fq), path(json), path(refDB)
    val params.database

    output:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"), path("${sample_name}.qc.xls"), path("map.depth.txt")
    tuple val(sample_name), path("${sample_name}.chromQC.xls"), path("${sample_name}.cov.png")
    tuple val(sample_name), path("${sample_name}.chromQC.html")

    when:
     refDB.size() > 0

    script:
    """
    ${params.python} ${params.script}/fasta_to_bed.py -i refDB.fasta -b refDB.bed
    ${params.bwa} index refDB.fasta

    ${params.bwa} mem -t $params.thread -R "@RG\\tID:lib_lane\\tLB:lib\\tSM:${sample_name}" -M \
        refDB.fasta ${sample_name}.clean.fq.gz | ${params.samtools} view -@ $params.thread -bS - |\
        ${params.samtools} sort -@ $params.thread - > ${sample_name}.sort.bam
    ${params.samtools} index -@ $params.thread ${sample_name}.sort.bam

    ${params.samtools} depth \
        -a -aa -b refDB.bed \
        ${sample_name}.sort.bam > map.depth.txt

    ${params.python} ${params.script}/quality_control.py \
        -bam ${sample_name}.sort.bam \
        -ref refDB.fasta \
        -json ${sample_name}.json \
        -out ${sample_name}.qc.xls \
        -depth map.depth.txt \
        -chrom ${sample_name}.chromQC.xls

    ${params.Rscript_cov} ${params.script}/genome_coverage.R map.depth.txt ${sample_name}.cov.png refDB.bed
    ${params.Rscript_html} ${params.script}/txt2html.R ${sample_name}.chromQC.xls
    """
}

// ==============================================
process cal_variation {
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.html", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.consensus.fasta", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.summaryQC.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.result.xls", mode: 'copy'
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.fastpQC.xls", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.sort.bam"), path("${sample_name}.sort.bam.bai"),path("${sample_name}.qc.xls"), path("map.depth.txt"), path(refDB),path("${sample_name}.json")

    output:
    tuple val(sample_name), path("${sample_name}.freeybayes.vcf"),path("${sample_name}_igv.html"), path("${sample_name}.consensus.fasta")
    tuple val(sample_name),path("${sample_name}.summaryQC.xls"), path("${sample_name}.result.xls"),path("${sample_name}.variantStat.xls"),path("${sample_name}.fastpQC.xls")

    when:
    refDB.size() > 0

    script:
    """
    ${params.samtools} faidx refDB.fasta

    ${params.freebayes} -C 5 -p 1 -F 0.8 -m 20 -q 20 --min-coverage 30 --limit-coverage 500 --pooled-continuous -f refDB.fasta ${sample_name}.sort.bam > ${sample_name}.freeybayes.raw.vcf

    ${params.perl} ${params.script}/decompose_mnp_v1.pl ${sample_name}.freeybayes.raw.vcf > ${sample_name}.freeybayes.vcf

    ${params.python} ${params.script}/low_region.py -i map.depth.txt  -d 1 -o ${sample_name}.low_cov.region

    ${params.bcftools} convert -Oz -o ${sample_name}.freeybayes.vcf.gz ${sample_name}.freeybayes.vcf
    ${params.bcftools} index -f ${sample_name}.freeybayes.vcf.gz
    ${params.bcftools} consensus -f refDB.fasta -o ${sample_name}.consensus.fasta.tmp -m ${sample_name}.low_cov.region -s ${sample_name} ${sample_name}.freeybayes.vcf.gz

    ${params.perl} ${params.script}/genome_stat.pl -vcf ${sample_name}.freeybayes.vcf -low_region ${sample_name}.low_cov.region -fa ${sample_name}.consensus.fasta.tmp -o ./ --sample $sample_name

    ${params.bedtools} genomecov -ibam ${sample_name}.sort.bam -bg > ${sample_name}.bedgraph

    ${params.python} ${params.script}/igv_html.py --fasta refDB.fasta --fai refDB.fasta.fai --vcf ${sample_name}.freeybayes.vcf --bedgraph ${sample_name}.bedgraph --out ${sample_name}_igv.html

    echo -e "seqName\tclade\tqc.overallStatus\tcoverage" > ${sample_name}.genotype.tsv

    for db in InfluenzaB/vic_ha InfluenzaB/yam_ha InfluenzaA/h3n2_ha InfluenzaA/h3n2_na InfluenzaA/h5n1 InfluenzaA/h1n1pdm_ha InfluenzaA/h1n1pdm_na;do
        ${params.nextclade} run -D ${params.database}/Nextclade/\$db --output-tsv \${db}.tsv ${sample_name}.consensus.fasta
        ${params.python} ${params.script}/pick_genotype.py -i \${db}.tsv -o \${db}.filter.tsv
        cat \${db}.filter.tsv >> ${sample_name}.genotype.tsv
    done

    ${params.python} ${params.script}/trans_result.py \
        -qc ${sample_name}.qc.xls \
        -json ${sample_name}.json \
        -genotype ${sample_name}.genotype.tsv \
        -vstat ${sample_name}.consensus_genome.stat \
        -outdir ./ \
        -sample $sample_name
    """
}

// ==============================================
process anno_mut {
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.html", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.freeybayes.vcf"), path(refDB)
    val params.genome

    output:
    tuple val(sample_name), path("${sample_name}.anno.html")

    when:
    refDB.size() > 0

    script:
    """
    ${params.python} ${params.script}/flu_build_snpeff.py \
        -fa refDB.fasta \
        -gff ${params.genome}/GFF/Influenza.gff \
        --map ${params.genome}/largeDB/segdedup.mapping.tsv \
        -outdir ./ \
        -java ${params.java} \
        -snpEff ${params.snpeff}

    ${params.java} -jar ${params.snpeff} ann -c snpEff.config -noStats -hgvs1LetterAa influenza ${sample_name}.freeybayes.vcf -canon > ${sample_name}.anno.vcf

    ${params.perl} ${params.script}/convert_snpeff_to_tab.pl ${sample_name}.anno.vcf > ${sample_name}.anno.xls

    ${params.Rscript_html} ${params.script}/txt2html.R ${sample_name}.anno.xls
    """
}

// ==============================================
process phylogenetic {
    publishDir "${params.outpath}/Result/${sample_name}", pattern: "*.png", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.qc.xls"), path(consensus)
    val params.database

    output:
    tuple val(sample_name), path("${sample_name}.seg4.phylip.png"), path("${sample_name}.seg6.phylip.png")

    when:
    consensus.size() > 0

    script:
    """
    ${params.python} ${params.script}/merge_seq_for_phlo.py \
        -i ${sample_name}.qc.xls \
        -c ${sample_name}.consensus.fasta \
        -d ${params.database}/Phylogenetic/PhyloDB.txt \
        -r ${params.genome}/largeDB/Influenza_wholeDB.fasta \
        -s ${params.database}/segment.txt \
        -m $sample_name \
        -o seg4.fasta \
        -n seg6.fasta

    ${params.mafft} --auto seg4.fasta > seg4.mafft.out
    ${params.fasttree} -nt -gtr seg4.mafft.out > ${sample_name}.seg4.fasttree

    ${params.mafft} --auto seg6.fasta > seg6.mafft.out
    ${params.fasttree} -nt -gtr seg6.mafft.out > ${sample_name}.seg6.fasttree

    ${params.Rscript_tree} ${params.script}/treeplot.R ${sample_name}.seg4.fasttree ${sample_name}.seg4.phylip.png $sample_name 4
    ${params.Rscript_tree} ${params.script}/treeplot.R ${sample_name}.seg6.fasttree ${sample_name}.seg6.phylip.png $sample_name 6
    """
}

// ==============================================
process phylogenetic_fa {
    publishDir "${params.outpath}/Result/", pattern: "*.png", mode: 'copy'

    input:
    val params.reffa
    val params.database

    output:
    tuple path("seg4.phylip.png"), path("seg6.phylip.png")

    script:
    """
    ${params.blastn} \
        -query $params.reffa -db ${params.genome}/largeDB/Influenza_wholeDB.fasta \
        -out blast.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qlen qstart qend sstrand slen sstart send evalue bitscore stitle" \
        -evalue 0.1 -max_target_seqs 15 -dust no -num_threads $params.thread

    ${params.python} ${params.script}/filter_blast_fa.py \
        -i blast.txt -o ./ -c $params.reffa \
        -r ${params.genome}/largeDB/Influenza_wholeDB.fasta \
        -g ${params.database}/GCF_seg_info.txt \
        -s ${params.database}/segment.txt \
        -d ${params.database}/Phylogenetic/PhyloDB.txt

    ${params.mafft} --auto seg4.fasta > seg4.mafft.out
    ${params.fasttree} -nt -gtr seg4.mafft.out > seg4.fasttree

    ${params.mafft} --auto seg6.fasta > seg6.mafft.out
    ${params.fasttree} -nt -gtr seg6.mafft.out > seg6.fasttree

    ${params.Rscript_tree} ${params.script}/treeplot.R seg4.fasttree seg4.phylip.png target 4
    ${params.Rscript_tree} ${params.script}/treeplot.R seg6.fasttree seg6.phylip.png target 6
    """
}

// ==============================================
process report1{
    publishDir "${params.outpath}/Result", pattern: "*.xls",mode: 'copy'

    input:
    path(chromQC)
    path(summaryQC)
    path(fastpQC)
    path(all_result)

    output:
    path "Report1.done"

    script:
    """
    ${params.python} ${params.script}/merge_qc.py -indir ${params.outpath}/Result -outdir ${params.outpath}/Result -tag 0
    touch Report1.done
    """
}

// ==============================================
process report{
    publishDir "${params.outpath}/Result",pattern: "*.xls" , mode: 'copy'

    input:
    path(chromQC)
    path(summaryQC)
    path(fastpQC)
    path(all_result)
    path(variation_stat)
    path(genome_stat)
    val(value)

    output:
    path "Report.done"

    script:
    """
    ${params.python} ${params.script}/merge_qc.py -indir ${params.outpath}/Result -outdir ${params.outpath}/Result
    mkdir FASTA

    ${params.python} ${params.script}/merge_all_seq_for_phlo.py \
        -i ${params.outpath}/Result/result.xls \
        -c ${params.outpath}/Result \
        -d ${params.database}/Phylogenetic/PhyloDB.txt \
        -r ${params.genome}/largeDB/Influenza_wholeDB.fasta \
        -s ${params.database}/segment.txt \
        -b ${params.script}/treeplot_all.R \
        -o FASTA \
        -rs ${params.Rscript_tree} \
        --mafft ${params.mafft} \
        --fasttree ${params.fasttree} 

    mkdir -p ${params.outpath}/Result/png
    cp FASTA/*png ${params.outpath}/Result/png

    ${params.python} ${params.script}/IF_auto_model_pmd.py --outdir ${params.outpath}/Result
    cp -r ${params.script}/IF_report/src ${params.outpath}/Result
    cp -r ${params.script}/IGV ${params.outpath}/Result

    touch Report.done
    """
}

// ==============================================
process json_result{
    publishDir "${params.outpath}", pattern: "*.json", mode: 'copy'
    input:
    val (value)
    val (value1)

    output:
    path "Result.json"

    script:
    """
    ${params.python} ${params.script}/json_result.py \
        -runid ${params.Runid} \
        -qc ${params.outpath}/Result/summaryQC.xls \
        -map ${params.outpath}/Result/result.xls \
        -stat ${params.outpath}/Result/variantStat.xls

    rm -rf ${params.outpath}/Result/report_files
    """
}
