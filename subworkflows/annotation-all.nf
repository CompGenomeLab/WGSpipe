include { FUNCOTATOR} from '../modules/funcotator.nf'
include { ANNOVAR } from '../modules/annovar.nf'
include { SNPEFF } from '../modules/snpeff.nf'

workflow ANNOTATION{

    take:
    tools
    ref
    fasta_index
    dict
    merge_vcf
    merge_vcf_idx
    snpeff_db
   
    main:
    
    if (tools.split(',').contains('funcotator')){
        FUNCOTATOR(ref,fasta_index, dict, merge_vcf,  merge_vcf_idx)
    }
    if (tools.split(',').contains('annovar')){
         ANNOVAR(ref,fasta_index, dict, merge_vcf,  merge_vcf_idx)
    }
    if (tools.split(',').contains('snpeff')){
         SNPEFF(ref,fasta_index, dict, merge_vcf,  merge_vcf_idx, snpeff_db)
    }
}