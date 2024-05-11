process FREEBAYES {
    tag "Executing Freebayes"
    publishDir "${params.outdir}/freebayes_call", mode: 'copy'


    input:
    path ref
    path applyed_bqsr_bam
    path bqsr_idx
    path fasta_index
    path dict

    output:
    path("freebayes_variants.vcf"), emit:fb //will produce a VCF file describing all SNPs, INDELs, and haplotype variants between the reference and aln.bam.
    
    script:
    """
    freebayes -f ${ref} ${applyed_bqsr_bam} > freebayes_variants.vcf
    """
}