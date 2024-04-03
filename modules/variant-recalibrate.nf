process VAR_RECAL {
    tag "Variant Recalibration"
    publishDir "${params.outdir}/haplotypecaller", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path variants
    path dbsnp
    path thousandG
    path dbsnp_idx
    path thousandG_idx
    path hapmap
    path hapmap_idx
    path omni
    path omni_idx

    output:
    path "${variants.baseName}_output.recal", emit: var_recal
    path "${variants.baseName}_output.tranches", emit: tranches
    path "${variants.baseName}_output.recal.idx", emit: var_recal_idx

    script:
    """
    gatk VariantRecalibrator \
            -R ${ref} \
            -V ${variants} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousandG} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode BOTH \
            -O ${variants.baseName}_output.recal \
            --tranches-file ${variants.baseName}_output.tranches
    """
}