process HAPPY{
    publishDir "${params.outdir}/happy_analysis/", mode: "copy"

    input:
    path ref
    path fasta_index //
    path refvcf //benchmark
    path refindex //benchmark index
    path refvcf_bed
    path vcf_file //my vcf

    output:
    tuple val("${vcf_file.getBaseName()}"), path("${vcf_file.getBaseName()}.*") 

    """
    /opt/hap.py/bin/hap.py \
    ${refvcf} ${vcf_file} \
    -f ${refvcf_bed} \
    -r ${ref} \
    --engine=vcfeval \
    --preprocess-truth \
    -o ${vcf_file.getBaseName()} \
    --pass-only \
    -l chr20
    """
}