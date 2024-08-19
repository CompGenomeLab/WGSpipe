process APPLY_VQSR {
    tag "Applying Variant Recalibration"
    publishDir "${params.outdir}/htvc_call", mode: 'copy'

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