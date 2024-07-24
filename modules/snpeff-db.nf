process SNPEFF_DB {

    input:
    val snpeff_db

    output:
    path ("\$SNPEFF_DATA_DIR"), emit:snpeff_cache
    
    script:
    """
    mkdir -p ${task.workDir}/snpeff-data

    export SNPEFF_DATA_DIR=${task.workDir}/snpeff-data

    snpEff download -dataDir \$SNPEFF_DATA_DIR ${snpeff_db}
    """
}