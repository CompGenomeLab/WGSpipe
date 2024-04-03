process APPLY_BQSR {
    tag "GATK applyBQSR underway"
    publishDir "${params.outdir}/apply-bqsr", mode: 'copy'

    input:
    path dedup_bam
    path ref
    path table
    path fasta_index
    path dict
    path interval

    output:
    path "*_recalibrated.bam", emit:applyed_bqsr_bam
    path "*_recalibrated.bam.bai", emit:bqsr_idx

    script:
    """
    gatk ApplyBQSR -I ${dedup_bam} \
    -R ${ref} \
    -L ${interval} \
    --bqsr-recal-file ${table} \
    -O ${dedup_bam.baseName}_recalibrated.bam
    
    gatk --java-options "-Xmx4g" BuildBamIndex -I ${dedup_bam.baseName}_recalibrated.bam -O ${dedup_bam.baseName}_recalibrated.bam.bai

    """
}