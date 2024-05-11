process BWA2_INDEX {

    input:
    path(ref)

    output:
    path("${refgenome}.0123"), emit: f0123
    path("${refgenome}.amb"), emit: amb
    path("${refgenome}.ann"), emit: ann
    path("${refgenome}.bwt.2bit.64"), emit: bwt_2bit_64
    path("${refgenome}.pac"), emit: pac

    script:
    """
    bwa-mem2 index ${ref}
    """
}

process BWA_MEM_2 { //need CPU and memory

    tag "BWA MEM 2 Alignment is underway"
    publishDir "${params.outdir}/${sample}-ALIGNED-2-bwa", mode: 'copy'


    input:
    path ref
    path bwa2_index
    tuple val(sample), path(reads)


    output:
    path "*.sam", emit: mem2_align

    script:
    """
    bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}" ${ref} ${reads[0]} ${reads[1]} > ${sample}_mem2_paired.sam
    """
}

