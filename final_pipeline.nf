params.genome = 'hg38' // Specify the genome name here

// Load parameters for the selected genome
params.ref = params.genomes[params.genome]?.ref
params.reads = params.genomes[params.genome]?.reads
params.site = params.genomes[params.genome]?.site
params.outdir = params.genomes[params.genome]?.outdir
params.dict = params.genomes[params.genome]?.dict
params.fasta_index = params.genomes[params.genome]?.fasta_index
params.site_idx = params.genomes[params.genome]?.site_idx
params.dbsnp = params.genomes[params.genome]?.dbsnp
params.thousandG = params.genomes[params.genome]?.thousandG
params.dbsnp_idx = params.genomes[params.genome]?.dbsnp_idx
params.thousandG_idx = params.genomes[params.genome]?.thousandG_idx
params.hapmap = params.genomes[params.genome]?.hapmap
params.hapmap_idx = params.genomes[params.genome]?.hapmap_idx
params.omni = params.genomes[params.genome]?.omni
params.omni_idx = params.genomes[params.genome]?.omni_idx

params.benchmark = "${launchDir}/new_data/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.bench_idx = "${launchDir}/new_data/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
params.bench_bed = "${launchDir}/new_data/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

log.info """\
      WGS PİPELİNE
================================
Reference        : ${params.ref}
Reads            : ${params.reads}
Output-folder    : ${params.outdir}/
Known-sites      : ${params.site}
"""

read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)

process fastqc {
  tag "Fastqc process underway"

  publishDir "${params.outdir}/fastqc_reports/", mode: 'copy'

  input:
  tuple val(sample), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit:fastqc_files

  script:
  """
  fastqc -t ${task.cpus} ${reads}
  """
}

process MULTIQC {
  tag "MULTIQC process underway"

  publishDir "${params.outdir}/multiqc/", mode: 'copy'
  
  input:
  tuple val(sample), path(reports)

  output:
  path("*.html")

  script:
  """
  multiqc ${reports}
  """
}

process BWA_INDEX {
    tag "Refrence Indexing by BWA is underway"
    publishDir "${params.outdir}/bwa-index/", mode: 'copy'

    input:
    path ref

    output:
    path("*") , emit: index

    script:
    """
    bwa index ${ref}
    """
}

process TRIM_FASTP{
    tag "Trimming with fastp process underway"

    publishDir "${params.outdir}/fastp_trimmed/", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_trimmed_1.fastq"), path("${sample}_trimmed_2.fastq"), emit: trimmed
    tuple val(sample), path('*.json'), emit: json
    tuple val(sample), path('*.html'), emit: html

    script:
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} \
    --thread ${task.cpus} \
    --detect_adapter_for_pe \
    --out1 ${sample}_trimmed_1.fastq --out2 ${sample}_trimmed_2.fastq \
    -j ${sample}_fastp.json -h ${sample}_fastp.html
    """
}

process BWA_ALIGNER { //need CPU and memory

    tag "BWA Alignment is underway"
    publishDir "${params.outdir}/${sample}-ALIGNED-bwa", mode: 'copy'


    input:
    path ref
    path index
    tuple val(sample), path(reads)

    output:
    path "*.sam", emit: aligned

    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}" ${ref} ${reads} > ${sample}_paired.sam
    """
}

process SAM_BAM {

    tag "Convert SAM to BAM"
    publishDir "${params.outdir}/sam-to-bam", mode: 'copy' // Output directory for BAM file

    input:
    path sam_file

    output:
    path "*_aligned.bam"
    path "*_sorted_aligned.bam", emit: bam

    script:
    """
    samtools view -bS ${sam_file} -o ${sam_file.baseName}_aligned.bam
    samtools sort -@ ${task.cpus} ${sam_file.baseName}_aligned.bam -o ${sam_file.baseName}_sorted_aligned.bam
    """
}

process MARK_DUP {

    tag "GatkMarkDuplicates underway"
    publishDir "${params.outdir}/mark_duplicated_aligned", mode: 'copy'

    input:
    path bam_file

    output:
    path "*_dedup.bam", emit: dedup_bam
    path('*bai') , emit: bai
    path('*sbi') , emit: sbi

    script:
    """
    gatk MarkDuplicatesSpark -I ${bam_file} -O ${bam_file.baseName}_dedup.bam --remove-sequencing-duplicates --conf 'spark.executor.cores=${task.cpus}'
    """
}

process BASE_RECAP {
    tag "GATK Base Recalibration underway"
    publishDir "${params.outdir}/recap_table", mode: 'copy'

    input:
    path dedup_bam
    path ref
    path fasta_index
    path dict
    path sites
    path sites_tbi

    output:
    path "*_recal_data.table", emit:table

    script:
    """
    gatk BaseRecalibratorSpark -I ${dedup_bam} -R ${ref}  --known-sites ${sites} -O ${dedup_bam.baseName}_recal_data.table
    """
}

process applyBQSR {
    tag "GATK applyBQSR underway"
    publishDir "${params.outdir}/apply_BQSR", mode: 'copy'

    input:
    path dedup_bam
    path ref
    path table
    path fasta_index
    path dict

    output:
    path "*_recalibrated.bam", emit:applyed_bqsr_bam
    path "*_recalibrated.bam.bai", emit:bqsr_idx

    script:
    """
    gatk ApplyBQSR -I ${dedup_bam} -R ${ref}  --bqsr-recal-file ${table} -O ${dedup_bam.baseName}_recalibrated.bam

    gatk --java-options "-Xmx4g" BuildBamIndex -I ${dedup_bam.baseName}_recalibrated.bam -O ${dedup_bam.baseName}_recalibrated.bam.bai
    """
}

process SUMMARY_METRICS {

    tag "Evaluating overall quality of the alignments"
    publishDir "${params.outdir}/SummeryMetrics", mode: 'copy'

    input:
    path ref
    path applyed_bqsr_bam

    output:
    path ('alignment_metrics.txt'), emit:Metrics
    path ('insert_size_metrics.txt'), emit:InsertSize
    path ('histogram.pdf'),  emit:histogram

    script:
    """
    gatk CollectAlignmentSummaryMetrics R=${ref} I=${applyed_bqsr_bam} O=alignment_metrics.txt
    gatk CollectInsertSizeMetrics INPUT=${applyed_bqsr_bam} OUTPUT=insert_size_metrics.txt \
    HISTOGRAM_FILE=histogram.pdf
    """

}

process HAPLOTYPECALLER {
    tag "Executing HaplotypeCaller"
    publishDir "${params.outdir}/variant_modules", mode: 'copy'

    input:
    path ref
    path applyed_bqsr_bam
    path fasta_index
    path dict

    output:
    path('htvc_variants.vcf'), emit:htvc

    script:
    """
    gatk HaplotypeCaller -R ${ref} -I ${applyed_bqsr_bam} -O htvc_variants.vcf
    """
}

process DEEPVARIANT {
    tag "Executing DeepVariant process"
    publishDir "${params.outdir}/variant_modules", mode: 'copy'

    input:
    path ref
    path applyed_bqsr_bam
    path applyed_bqsr_bai
    path fasta_index
    path dict

    output:
    path('dv_variants.vcf'), emit: dv

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${ref} \
        --reads=${applyed_bqsr_bam} \
        --output_vcf=dv_variants.vcf \
        --num_shards=${task.cpus} \
        --regions chr20 \

    """
}

process VAR_RECAL {
    tag "Variant Recalibration"
    publishDir "${params.outdir}/haplotypecaller", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path variants
    path dbsnp
    path thousandG
    path dbsnp_idx
    path thousandG_idx
    path hapmap
    path hapmap_idx
    path omni
    path omni_idx

    output:
    path "${variants.baseName}_output.recal", emit: var_recal
    path "${variants.baseName}_output.tranches", emit: tranches
    path "${variants.baseName}_output.recal.idx", emit: var_recal_idx

    script:
    """
    gatk VariantRecalibrator \
            -R ${ref} \
            -V ${variants} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousandG} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode BOTH \
            -O ${variants.baseName}_output.recal \
            --tranches-file ${variants.baseName}_output.tranches
    """
}

process apply_VQSR {
    tag "Applying Variant Recalibration"
    publishDir "${params.outdir}/haplotypecaller", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path variants
    path var_recal
    path tranches
    path var_recal_idx

    output:
    path "${variants.baseName}_recalibrated_variant.vcf", emit:htvc_recalibrated
    path "*.idx", emit:htvc_index_recalibrated

    script:
    """
    gatk ApplyVQSR \
        -R ${ref} \
        -V ${variants}\
        --recal-file ${var_recal} \
        --tranches-file ${tranches} \
        -mode BOTH \
        -ts-filter-level 99.0 \
        -O ${variants.baseName}_recalibrated_variant.vcf
    """
}

process VARIANT_FILTER {

    publishDir "${params.outdir}/haplotypecaller", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path recal_variants
    path recal_idx_variants

    output:
    path "htvc_variant_filtered.vcf", emit:htvc_filtered
    path "*.idx"


    script:
    """
    gatk VariantFiltration \
            -R ${ref} \
            -V ${recal_variants} \
            -O htvc_variant_filtered.vcf \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 60.0" \
            -filter-name "MQ_filter" -filter "MQ < 40.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    """
}

process SNP {
    tag "Selecting SNP variants"

    publishDir "${params.outdir}/haplotypecaller", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path variants

    output:
    path "htvc_snp_filtered_variants.vcf", emit: htvc_snp

    script:
    """
    gatk SelectVariants -R ${ref} -V ${variants} --select-type-to-include SNP -O htvc_snp_filtered_variants.vcf
    """
}

//selecting INDEL from filtered vcf file
process INDEL {
    tag "Selecting INDEL variants"

    publishDir "${params.outdir}/haplotypecaller", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path variants

    output:
    path "htvc_indel_filtered_variants.vcf", emit:htvc_indel

    script:
    """
    gatk SelectVariants -R ${ref} -V ${variants} --select-type-to-include INDEL -O htvc_indel_filtered_variants.vcf
    """
}

process SNP_DV {
    tag "Selecting SNP variants"

    publishDir "${params.outdir}/deepvariant", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path variants

    output:
    path "dv_snp_filtered_variants.vcf", emit: dv_snp

    script:
    """
    gatk SelectVariants -R ${ref} -V ${variants} --select-type-to-include SNP -O dv_snp_filtered_variants.vcf
    """
}

process INDEL_DV {
    tag "Selecting INDEL variants"

    publishDir "${params.outdir}/deepvariant", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path variants

    output:
    path "dv_indel_filtered_variants.vcf", emit:dv_indel

    script:
    """
    gatk SelectVariants -R ${ref} -V ${variants} --select-type-to-include INDEL -O dv_indel_filtered_variants.vcf
    """
}

process MERGE_VCFS_HTVC {

    tag "Merge HTVC Variant Files"

    publishDir "${params.outdir}/htvc_merge_variants", mode: 'copy'

    input:
    path snp
    path indel

    output:
    path "htvc_merged_variants.vcf", emit:htvc_merged
    path "*.vcf.idx", emit:htvc_merged_idx

    script:
    """
    gatk MergeVcfs -I ${snp} -I ${indel} -O htvc_merged_variants.vcf
    gatk IndexFeatureFile -I htvc_merged_variants.vcf
    """
}

process MERGE_VCFS_DV {

    tag "Merge DV Variant Files"

    publishDir "${params.outdir}/dv_merge_variants", mode: 'copy'

    input:
    path snp
    path indel

    output:
    path "dv_merged_variants.vcf", emit:dv_merged
    path "*.vcf.idx", emit:dv_merged_idx

    script:
    """
    gatk MergeVcfs -I ${snp} -I ${indel} -O dv_merged_variants.vcf
    gatk IndexFeatureFile -I dv_merged_variants.vcf
    """
}

process FUNCOTATOR_ANNOTATION_HTVC {

    tag "Selecting INDEL variants"

    publishDir "${params.outdir}/htvc_funcotator_annotation", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx

    output:
    path "htvc_annotation.vcf"

    script:
    """
    gatk Funcotator \
    -R ${ref} \
    -V ${merge_variant}\
    -O htvc_annotation.vcf \
    --output-file-format VCF \
    --data-sources-path ${launchDir}/funcotator_data/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38
    """
}

process FUNCOTATOR_ANNOTATION_DV {

    tag "Funcotator annotation on dv"

    publishDir "${params.outdir}/dv_funcotator_annotation", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx

    output:
    path "dv_funcotator_annotation.vcf"

    script:
    """
    gatk Funcotator \
    -R ${ref} \
    -V ${merge_variant}\
    -O dv_funcotator_annotation.vcf \
    --output-file-format VCF \
    --data-sources-path ${launchDir}/funcotator_data/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38
    """
}

process ANNOVAR_HTVC{
    tag "Annotating htvc variants with ANNOVAR"
    publishDir "${params.outdir}/annovar_htvc", mode: 'copy'

    input:
    path ref
    path fasta_index // Index of ref
    path dict // Dictionary of ref
    path htvc_variant // Input variant file (VCF format)
    path htvc_variant_idx

    output:
    path "*.hg38_multianno.vcf"
    path "*.hg38_multianno.txt"
    path "*.avinput"

    script:
    """
    table_annovar.pl ${htvc_variant} ${launchDir}/humandb/ \
        --buildver hg38 \
        --out ${htvc_variant.baseName}_annovar \
        --remove \
        --protocol refGene,gnomad_genome,gnomad_exome,avsnp150,clinvar_20221231 \
        --operation g,f,f,f,f \
        --nastring . \
        --vcfinput \
        --thread ${task.cpus}
    """
}

process ANNOVAR_DV {
    tag "Annotating dv variants with ANNOVAR"
    publishDir "${params.outdir}/annovar_dv", mode: 'copy'

    input:
    path ref
    path fasta_index // Index of ref
    path dict // Dictionary of ref
    path dv_variant // Input variant file (VCF format)
    path dv_variant_idx

    output:
    path "*.hg38_multianno.vcf"
    path "*.hg38_multianno.txt"
    path "*.avinput"

    script:
    """
    table_annovar.pl ${dv_variant} ${launchDir}/humandb/ \
        --buildver hg38 \
        --out ${dv_variant.baseName}_annovar \
        --remove \
        --protocol refGene,gnomad_genome,gnomad_exome,avsnp150,clinvar_20221231 \
        --operation g,f,f,f,f \
        --nastring . \
        --vcfinput \
        --thread ${task.cpus}
    """
}

process process_happy{
    publishDir "${params.outdir}/compare_w_happy/", mode: "copy"

    tag {vcf_file.getBaseName()}

    input:
    path ref
    path fasta_index //
    path refvcf //benchmark
    path refindex //benchmark index
    path refvcf_bed
    path vcf_file //my vcf

    output:
    tuple val("${vcf_file.BaseName()}"), path("${vcf_file.BaseName()}.*") 

    """
    /opt/hap.py/bin/hap.py \
    ${refvcf} ${vcf_file} \
    -f ${refvcf_bed} \
    -r ${ref} \
    --engine=vcfeval \
    --preprocess-truth \
    -o ${vcf_file.BaseName()}
    --pass-only \
    -l chr20
    """
}

workflow pipeline {

    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    fastqc.out.fastqc_files.view()
    MULTIQC(fastqc.out.fastqc_files)
    
    BWA_INDEX(params.ref)
    TRIM_FASTP(read_pairs_ch)
    BWA_ALIGNER(params.ref, BWA_INDEX.out.index, TRIM_FASTP.out.trimmed)
    SAM_BAM(BWA_ALIGNER.out.aligned)

    MARK_DUP(SAM_BAM.out.bam)
    BASE_RECAP(MARK_DUP.out.dedup_bam, params.ref, params.fasta_index, params.dict, params.site, params.site_idx)
    applyBQSR(MARK_DUP.out.dedup_bam, params.ref, BASE_RECAP.out.table, params.fasta_index, params.dict)
    SUMMARY_METRICS(params.ref, applyBQSR.out.applyed_bqsr_bam)

    HAPLOTYPECALLER(params.ref, applyBQSR.out.applyed_bqsr_bam,params.fasta_index, params.dict)
    DEEPVARIANT(params.ref, applyBQSR.out.applyed_bqsr_bam, applyBQSR.out.bqsr_idx, params.fasta_index, params.dict)

    VAR_RECAL(params.ref, params.fasta_index, params.dict, HAPLOTYPECALLER.out.htvc, params.dbsnp, params.thousandG, params.dbsnp_idx, params.thousandG_idx, params.hapmap, params.hapmap_idx, params.omni, params.omni_idx)

    apply_VQSR(params.ref,params.fasta_index, params.dict, HAPLOTYPECALLER.out.htvc, VAR_RECAL.out.var_recal, VAR_RECAL.out.tranches, VAR_RECAL.out.var_recal_idx)

    VARIANT_FILTER(params.ref, params.fasta_index, params.dict, apply_VQSR.out.htvc_recalibrated, apply_VQSR.out.htvc_index_recalibrated)

    SNP(params.ref, params.fasta_index, params.dict, VARIANT_FILTER.out.htvc_filtered)
    INDEL(params.ref, params.fasta_index, params.dict, VARIANT_FILTER.out.htvc_filtered)

    SNP_DV(params.ref, params.fasta_index, params.dict, DEEPVARIANT.out.dv)
    INDEL_DV(params.ref, params.fasta_index, params.dict, DEEPVARIANT.out.dv)

    MERGE_VCFS_HTVC(SNP.out.htvc_snp, INDEL.out.htvc_indel)
    MERGE_VCFS_DV(SNP_DV.out.dv_snp, INDEL_DV.out.dv_indel)

    FUNCOTATOR_ANNOTATION_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    FUNCOTATOR_ANNOTATION_DV(params.ref, params.fasta_index, params.dict, MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)

    ANNOVAR_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    ANNOVAR_DV(params.ref, params.fasta_index, params.dict,  MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)

    process_happy(params.ref, params.fasta_index, params.benchmark, params.bench_idx, params.bench_bed, DEEPVARIANT.out.dv)

}

workflow case_study{

    params.bam = "/cta/users/baharsevgin/test/input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"
    params.index = "/cta/users/baharsevgin/test/input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai"
    params.outdir = "${launchDir}/case-study" 

    
    DEEPVARIANT(params.ref, params.bam, params.index, params.fasta_index, params.dict)
    process_happy(params.ref, params.fasta_index, params.benchmark, params.bench_idx, params.bench_bed, DEEPVARIANT.out.dv)

}

workflow {
    pipeline()
    case_study()
}


