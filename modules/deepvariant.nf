process DEEPVARIANT {
    tag "Executing DeepVariant process"
    publishDir "${params.outdir}/deepvariant_call", mode: 'copy'

    input:
    path ref
    path applyed_bqsr_bam
    path applyed_bqsr_bai
    path fasta_index
    path dict
    path wes_bed

    output:
    path('dv_variants.vcf'), emit: dv

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${ref} \
        --reads=${applyed_bqsr_bam} \
        --output_vcf=dv_variants.vcf \
        --regions ${wes_b} \
        --num_shards=${task.cpus}
    """
}
