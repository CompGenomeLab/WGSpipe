process SNPEFF {
    tag "SnpEFF annotation"

    publishDir "${params.outdir}/annotation/${merge_variant.baseName}_snpeff", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx
    path snpeff_db

    output:
    path "*.ann.vcf"
    path "*.csv"

    script:
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        -v ${snpeff_db} \\
        -dataDir \${PWD}/tmp \\
        -csvStats ${merge_variant.baseName}_variant_snpeff.csv\\
        ${merge_variant} \\
        > ${merge_variant.baseName}_variant_snpeff.ann.vcf 
    """
}

// snpEff -v -download hg38  -i ${merge_variant} -o vcf -csvStats ${merge_variant.baseName}_variant_snpeff.csv > ${merge_variant.baseName}_variant_snpeff.ann.vcf

