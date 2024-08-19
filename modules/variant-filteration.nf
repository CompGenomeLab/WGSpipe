process VARIANT_FILTER {

    publishDir "${params.outdir}/htvc_call", mode: 'copy'

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