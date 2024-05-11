process run_mutect2 {
    input:
    path bam // Input BAM file
    path reference // Reference genome
    path vcf // VCF for known variants
   

    output:
    path(*.vcf) // Mutect2 output VCF

    script:
    """
    gatk Mutect2 \
        -R ${reference} \
        -I ${bam} \
        --germline-resource ${vcf} \
        -O ${bam.baseName}_mutect_variant.vcf
    """
}

