process ENSEMBL_VEP {
    publishDir "${params.outdir}/annotation/ensembl_vep/${vcf.baseName}_vep", mode: 'copy'

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
    path("*.vcf"), emit: vcf
    path("*.summary.html"), optional: true, emit: report

    script:
    """
    vep \\
	-i $vcf \\
        -o  ${vcf.baseName}_ensembl_vep.vcf \\
        --stats_file  ${vcf.baseName}_ensembl_vep.summary.html \\
        --fasta ${ref} \\
        --assembly ${genome} \\
        --species ${species} \\
        --cache \\
        --cache_version ${cache_version} \\
        --dir_cache ${launchDir}/vep-data \\
        --fork ${task.cpus}
    """
}
