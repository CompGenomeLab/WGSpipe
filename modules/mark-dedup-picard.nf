process MARK_DEDUP_PICARD {

    tag "GatkMarkDuplicates underway"
    publishDir "${params.outdir}/mark_duplicates_picard_gatk", mode: 'copy'

    input:
    path bam_file
    path bam_file_idx    

    output:
    path ("*_dedup.bam"), emit: dedup_bam
    path('*bai') , emit: bai
    path('*metrics')
    
    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates \
      -I ${bam_file} \
      -O ${bam_file.baseName}_dedup.bam \
      -M ${bam_file.baseName}_dedup.metrics \
      --TMP_DIR . \
      --CREATE_INDEX
    """
}
