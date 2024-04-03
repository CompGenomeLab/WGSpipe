process TRIM_FASTP{
    tag "Trimming with fastp process underway"

    publishDir "${params.outdir}/fastp_trimmed/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    
    
    output:
    tuple val(sample), path("${sample}_trimmed_1.fastq"), path("${sample}_trimmed_2.fastq"), emit: trimmed
    tuple val(sample), path('*.json'), emit: json
    tuple val(sample), path('*.html'), emit: html

    script:
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} \
    --detect_adapter_for_pe \
    --thread ${task.cpus} \
    --out1 ${sample}_trimmed_1.fastq --out2 ${sample}_trimmed_2.fastq \
    -j ${sample}_fastp.json -h ${sample}_fastp.html
    """
}