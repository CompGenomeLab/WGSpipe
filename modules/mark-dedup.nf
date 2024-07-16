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
    def avail_mem = (task.memory.mega*0.8).intValue()
    """
    gatk MarkDuplicatesSpark \
    -I ${bam_file} \
    -O ${bam_file.baseName}_dedup.bam \
    -M ${bam_file.baseName}_dedup.metrics \
    --conf 'spark.executor.cores=${task.cpus}'
    """
}

 gatk MarkDuplicatesSpark \
    -I HG001.novaseq.pcr-free.30x.dedup.grch38.bam \
    -O test_dedup.bam \
    -M test_dedup.metrics \
    --tmp-dir /cta/users/baharsevgin/tmp
    --conf 'spark.executor.cores=6'
    

/cta/users/baharsevgin/ENS491-pipeline/results_HG003_WGS/sam-to-bam