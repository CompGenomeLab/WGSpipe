process BWA_MEM { //need CPU and memory

    tag "BWA Alignment is underway"
    publishDir "${params.outdir}/${sample}-ALIGNED-bwa", mode: 'copy'


    input:
    path ref
    path index
    tuple val(sample), path(reads)

    output:
    path "*.sam", emit: aligned

    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}" ${ref} ${reads[0]} ${reads[1]} > ${sample}_paired.sam
    """
}