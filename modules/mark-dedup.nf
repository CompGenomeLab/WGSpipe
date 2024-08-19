process MARK_DEDUP {

    tag "GatkMarkDuplicates underway"
    publishDir "${params.outdir}/mark_duplicates_gatk", mode: 'copy'

    input:
    path bam_file
    path bam_file_idx    

    output:
    path ("*_dedup.bam"), emit: dedup_bam
    path('*bai') , emit: bai
    path('*sbi') , emit: sbi
    path('*metrics')
    
    script:
    """
    gatk MarkDuplicatesSpark \
    -I ${bam_file} \
    -O ${bam_file.baseName}_dedup.bam \
    -M ${bam_file.baseName}_dedup.metrics \
    --conf 'spark.executor.cores=${task.cpus}'
    """
}

 