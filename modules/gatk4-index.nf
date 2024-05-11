process GATK_INDEX{
    tag "Indexing VCF files" 
    publishDir "${params.outdir}/gatk4-index-vcffiles", mode: 'copy'

    input:
    path vcf

    output:
    path("*.vcf.idx"), emit:vcf_idx

    script:
    """
    gatk IndexFeatureFile -I ${vcf}
    """
}