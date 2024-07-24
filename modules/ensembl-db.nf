process ENSEMBLVEP_DOWNLOAD {

    input:
    val (assembly) 
    val (species)
    val (cache_version)

    output:
    tuple val(meta), path(prefix), emit: cache
    path "versions.yml"          , emit: versions

    script:
    """
    vep_install \\
        --CACHEDIR $prefix \\
        --SPECIES $species \\
        --ASSEMBLY $assembly \\
        --CACHE_VERSION $cache_version \\
        $args

    """
}