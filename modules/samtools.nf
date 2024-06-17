process SAM_CONVERTER {

    tag "Convert SAM to BAM"
    publishDir "${params.outdir}/sam-to-bam", mode: 'copy' // Output directory for BAM file

    input:
    path sam_file

    output:
    path "*_aligned.bam"
    path "*_sorted_aligned.bam", emit: bam
    path "*bai", emit: bam_bai
    path "*flagstat"

    script:
    """
    samtools view -bS ${sam_file} -o ${sam_file.baseName}_aligned.bam
    samtools sort -@ ${task.cpus} ${sam_file.baseName}_aligned.bam -o ${sam_file.baseName}_sorted_aligned.bam
    samtools index ${sam_file.baseName}_sorted_aligned.bam
    samtools flagstat ${sam_file.baseName}_sorted_aligned.bam >	${sam_file.baseName}.flagstat
    """
}