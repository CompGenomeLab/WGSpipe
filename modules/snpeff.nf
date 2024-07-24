process SNPEFF {
    tag "SnpEFF annotation"

    publishDir "${params.outdir}/annotation/${merge_variant.baseName}_snpeff", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx
    val snpeff_db
    

    output:
    path "*.ann.vcf"
    path "*.csv"

    script:
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
	avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    snpEff \\
	-Xmx${avail_mem}M \\
        -v ${snpeff_db} \\
        -dataDir ${launchDir}/snpeff-data \\
        -csvStats ${merge_variant.baseName}_variant_snpeff.csv\\
        ${merge_variant} \\
        > ${merge_variant.baseName}_variant_snpeff.ann.vcf
    """
}

// snpEff -v -download hg38  -i ${merge_variant} -o vcf -csvStats ${merge_variant.baseName}_variant_snpeff.csv > ${merge_variant.baseName}_variant_snpeff.ann.vcf

