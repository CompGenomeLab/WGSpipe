process DEEPVARIANT {
    tag "Executing DeepVariant process"
    publishDir "${params.outdir}/deepvariant_call", mode: 'copy'

    input:
    path ref
    path applyed_bqsr_bam
    path applyed_bqsr_bai
    path fasta_index
    path dict

    output:
    path('dv_variants.vcf'), emit: dv

    script:
    def modelType = params.wes ? "WES" : "WGS"
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=${modelType} \
        --ref=${ref} \
        --reads=${applyed_bqsr_bam} \
        --output_vcf=dv_variants.vcf \
        --regions ${launchDir}/data/idt_capture.grch38.bed \
        --num_shards=${task.cpus}
    """
}
