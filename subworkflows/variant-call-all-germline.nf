include { HAPLOTYPECALLER } from '../modules/haplotypecaller.nf'
include { DEEPVARIANT } from '../modules/deepvariant.nf'
include { FREEBAYES } from '../modules/freebayes.nf'

workflow VARIANT_CALLING{
    take:
    tools
    ref
    fasta_index
    dict
    applyed_bqsr_bam
    applyed_bqsr_bam_idx
    dbsnp
    dbsnp_idx
    interval

    main:

    vcf_htvc = Channel.empty()
    vcf_dv = Channel.empty()
    vcf_fb = Channel.empty()

    if(tools.split(',').contains('haplotypecaller')){
        HAPLOTYPECALLER(ref, applyed_bqsr_bam, applyed_bqsr_bam_idx, fasta_index, dict, dbsnp, dbsnp_idx, interval)
        vcf_htvc = HAPLOTYPECALLER.out.htvc
    }

    if(tools.split(',').contains('deepvariant')){
        DEEPVARIANT(ref, applyed_bqsr_bam, applyed_bqsr_bam_idx, fasta_index, dict)
        vcf_dv = DEEPVARIANT.out.dv
    }
    
    if(tools.split(',').contains('freebayes')){
        FREEBAYES(ref, applyed_bqsr_bam, applyed_bqsr_bam_idx, fasta_index, dict)
        vcf_fb = FREEBAYES.out.fb
    }
    
    emit:
    vcf_htvc
    vcf_dv
    vcf_fb
}

