include { HAPLOTYPECALLER } from '../modules/haplotypecaller.nf'
include { DEEPVARIANT } from '../modules/deepvariant.nf'
include { FREEBAYES } from '../modules/freebayes.nf'

workflow VARIANT_CALLING{
    take:

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
    vcf_all = Channel.empty()

    HAPLOTYPECALLER(ref, applyed_bqsr_bam, applyed_bqsr_bam_idx, fasta_index, dict, dbsnp, dbsnp_idx, interval)
    vcf_htvc = HAPLOTYPECALLER.out.htvc

    DEEPVARIANT(ref, applyed_bqsr_bam, applyed_bqsr_bam_idx, fasta_index, dict)
    vcf_dv = DEEPVARIANT.out.dv

    FREEBAYES(ref, applyed_bqsr_bam, applyed_bqsr_bam_idx, fasta_index, dict)
    vcf_fb = FREEBAYES.out.fb
 

    vcf_all = Channel.empty().mix(
        vcf_htvc,
        vcf_dv,
        vcf_fb
    )

    emit:
    vcf_htvc
    vcf_dv
    vcf_fb
    vcf_all
}

process run_manta {
    input:
    file bam // Input BAM file
    file reference // Reference genome
    string sample_name // Sample name for output naming
    
    output:
    file 'results/manta/*' // Manta results folder

    script:
    """
    configManta.py \
        --bam ${bam} \
        --referenceFasta ${reference} \
        --runDir results/manta_${sample_name}
        
    results/manta_${sample_name}/runWorkflow.py -m local -j 4 # Adjust job count based on your system
    """
}

// Define process for Mutect2
process run_mutect2 {
    input:
    file bam // Input BAM file
    file reference // Reference genome
    file vcf // VCF for known variants
    string output_vcf // Output VCF file name

    output:
    file output_vcf // Mutect2 output VCF

    script:
    """
    gatk Mutect2 \
        -R ${reference} \
        -I ${bam} \
        --germline-resource ${vcf} \
        -O ${output_vcf}
    """
}

// Define process for CNVkit
process run_cnvkit {
    input:
    file bam // Input BAM file
    file reference // Reference genome
    string output_dir // Output directory name

    output:
    file "${output_dir}/*" // CNVkit results folder

    script:
    """
    cnvkit.py batch ${bam} \
        --output-dir ${output_dir} \
        --reference ${reference} \
        --scatter \
        --diagram
    """
}