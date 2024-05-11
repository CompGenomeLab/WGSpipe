process SNPEFF {
    tag "SnpEFF annotation"

    publishDir "${params.outdir}/snpeff_annotation", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx

    output:
    path "*.ann.vcf"
    path "*.csv"

    script:
    """
    snpEff -v -download hg38  -i ${merge_variant} -o vcf -csvStats ${merge_variant.baseName}_variant_snpeff.csv > ${merge_variant.baseName}_variant_snpeff.ann.vcf
    """
}

