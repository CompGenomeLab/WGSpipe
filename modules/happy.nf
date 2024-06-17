process HAPPY{
    publishDir "${params.outdir}/happy_analysis/", mode: "copy"
    
    tag {vcf_file.getBaseName()}

    input:
    path ref 
    path fasta_index //
    path refvcf //benchmark
    path refindex //benchmark index
    path refvcf_bed
    path vcf_file 
    path wes_bed

    output:
    tuple val("${vcf_file.getBaseName()}"), path("${vcf_file.getBaseName()}.*") 

    """
    /opt/hap.py/bin/hap.py \
    ${refvcf} \
    ${vcf_file} \
    -r ${ref} \
    -f ${refvcf_bed} \
    -T ${wes_bed} \
    -o ${vcf_file.getBaseName()} \
    --engine=vcfeval \
    --pass-only \
    """
}
