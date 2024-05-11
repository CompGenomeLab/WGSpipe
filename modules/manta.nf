process manta {
    tag "Variant call with MANTA" 
    publishDir "${params.outdir}/variant-call-manta", mode: 'copy'

    input:
    path bam // Input BAM file
    path bam_bai
    path reference // Reference genome

    output:
    
    
    script:
    """
    configManta.py  --bam ${bam} --referenceFasta ${reference} --runDir .
    ./runWorkflow.py -m local -j ${task.cpus}
    """
}