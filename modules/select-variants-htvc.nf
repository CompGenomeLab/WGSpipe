process SELECT_SNP_HTVC {
    tag "Selecting SNP variants"

    publishDir "${params.outdir}/select-variant-htvc", mode: 'copy'

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
process SELECT_INDEL_HTVC {
    tag "Selecting INDEL variants"

    publishDir "${params.outdir}/select-variant-htvc", mode: 'copy'

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