include { FUNCOTATOR} from '../modules/funcotator.nf'
include { ANNOVAR } from '../modules/annovar.nf'
include { SNPEFF } from '../modules/snpeff.nf'

workflow ANNOTATION{

    take:
    ref
    fasta_index
    dict
    merge_vcf
    merge_vcf_idx
   
    main:
    
    FUNCOTATOR(ref,fasta_index, dict, merge_vcf,  merge_vcf_idx)
    ANNOVAR(ref,fasta_index, dict, merge_vcf,  merge_vcf_idx)
    SNPEFF(ref,fasta_index, dict, merge_vcf,  merge_vcf_idx)

}