process ENSEMBL_VEP {
    publishDir "${params.outdir}/annotation/${merge_variant.baseName}_ensembl_vep", mode: 'copy'

    input:
    path  vcf //generated variant files
    val   genome //GRCh38
    val   species //homo_sapiens
    val   cache_version //GRCh38 version
    path  ensembl_cache //directory of the database
    path  ref
    
    output:
    path("*.vcf.gz"), emit: vcf
    path("*.tab.gz"), emit: tab
    path("*.json.gz"), emit: json
    path "*.summary.html", emit: report

    script:
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    """
    vep \\
        -i $vcf \\
        -o  ${vcf.baseName}_ensembl_vep.${file_extension}.gz \\
        --fasta ${ref} \\
        --assembly ${genome} \\
        --species ${species} \\
        --cache \\
        --cache_version ${cache_version} \\
        --dir_cache ${ensembl_cache} \\
        --fork ${task.cpus}
    """
}