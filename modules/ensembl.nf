process ENSEMBL_VEP {
    publishDir "${params.outdir}/annotation/${vcf.baseName}_ensembl_vep", mode: 'copy'

    input:
    path  vcf //generated variant files
    path  vcf_idx
    val   genome //GRCh38
    val   species //homo_sapiens
    val   cache_version //GRCh38 version
    path  ref
    path  fasta_index
    path  dict
    
    output:
    path("*.vcf.gz"), emit: vcf
    path ("*.summary.html"), optional: true, emit: report

    script:
    """
    vep \\
        -i $vcf \\
        --stats_file  ${vcf.baseName}_ensembl_vep.summary.html \\
        -o ${vcf.baseName}_ensembl_vep.vcf.gz \\
        --fasta ${ref} \\
        --database \\
        --assembly ${genome} \\
        --species ${species} \\
        --cache_version ${cache_version} \\
        --fork ${task.cpus}
    """
}