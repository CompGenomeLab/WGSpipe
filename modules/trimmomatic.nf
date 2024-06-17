process TRIMMOMATIC{
    tag "Indexing VCF files" 
    publishDir "${params.outdir}/trimmomatic", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    val trimmomatic_adapters 
    val trimmomatic_adapters_param 
    val trimmomatic_window_len 
    val trimmomatic_window_val 
    val trimmomatic_min_len 

    output:
    tuple val(sample), path("${sample}_trimmed_{R1,R2}.fq.gz"), emit:trimmomatic_fastq
    tuple val(sample), path("${sample}_unpaired_{R1,R2}.fq.gz"), emit:unpaired

    script:
    """
     trimmomatic PE -phred33 -threads ${task.cpus} -summary ${reads}  "${sample}_trimmed_R1.fq.gz" "${sample}_unpaired_R1.fq.gz" \
     "${sample}_trimmed_R2.fq.gz" "${sample}_unpaired_R2.fq.gz" \
     ILLUMINACLIP:${params.trimmomatic_adapters}:${params.trimmomatic_adapters_param} \
     SLIDINGWINDOW:${params.trimmomatic_window_len}:${params.trimmomatic_window_val} \
     MINLEN:${params.trimmomatic_min_len} 2> ${name}.log
    """

}