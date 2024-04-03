process SELECT_SNP_DV {
    tag "Selecting SNP variants"

    publishDir "${params.outdir}/select-variant-dv", mode: 'copy'

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

//selecting INDEL from filtered vcf file
process SELECT_INDEL_DV {
    tag "Selecting INDEL variants"

    publishDir "${params.outdir}/select-variant-dv", mode: 'copy'

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