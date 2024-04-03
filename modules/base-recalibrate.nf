process BASE_RECAP {
    tag "GATK Base Recalibration underway"
    publishDir "${params.outdir}/base-recalibrator", mode: 'copy'

    input:
    path dedup_bam
    path dedup_bam_idx
    path ref
    path fasta_index
    path dict
    path sites
    path sites_tbi
    path interval
    path dbsnp
    path dbsnp_idx
    path indel
    path indel_idx

    output:
    path "*_recal_data.table", emit:table

    script:
    """
    gatk BaseRecalibratorSpark -I ${dedup_bam} -R ${ref} \
     --known-sites ${sites} \
     --known-sites ${dbsnp} \
     --known-sites ${indel}  \
     --intervals ${interval} \
     -O ${dedup_bam.baseName}_recal_data.table
    """
}