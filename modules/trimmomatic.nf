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
    tuple val(sample_id), path("trimmed_*_paired.fq.gz"), path("trimmed_*_unpaired.fq.gz") 

    script:
    """
     trimmomatic PE -phred33 ${reads} -threads ${task.cpus} "${sample}_filtered_R1.fastq" "${sample}_unpaired_R1.fastq" \
     "${sample}_filtered_R2.fastq" "${sample}_unpaired_R2.fastq" \
     ILLUMINACLIP:${params.trimmomatic_adapters}:${params.trimmomatic_adapters_param} \
     SLIDINGWINDOW:${params.trimmomatic_window_len}:${params.trimmomatic_window_val} \
     MINLEN:${params.trimmomatic_min_len} 2> ${name}.log
    
    """


}