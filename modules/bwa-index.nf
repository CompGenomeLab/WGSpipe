process BWA_INDEX {
    tag "Refrence Indexing by BWA is underway"
    publishDir "${params.outdir}/bwa-index/", mode: 'copy'

    input:
    path ref

    output:
    path("*") , emit: index

    script:
    """
    bwa index ${ref}
    """
}