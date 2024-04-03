process HAPLOTYPECALLER {
    tag "Executing HaplotypeCaller"
    publishDir "${params.outdir}/htvc_call", mode: 'copy'

    input:
    path ref
    path applyed_bqsr_bam
    path bqsr_idx
    path fasta_index
    path dict
    path dbsnp
    path dbsnp_idx
    path interval

    output:
    path('htvc_variants.vcf'), emit:htvc

    script:
    """
    gatk HaplotypeCaller -R ${ref} \
    -I ${applyed_bqsr_bam} \
    -D ${dbsnp} \
    -L ${interval} \
    -O htvc_variants.vcf
    """
}