process FUNCOTATOR_ANNOTATION_HTVC {

    tag "Selecting INDEL variants"

    publishDir "${params.outdir}/htvc_funcotator_annotation", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx

    output:
    path "htvc_annotation.vcf"

    script:
    """
    gatk Funcotator \
    -R ${ref} \
    -V ${merge_variant}\
    -O htvc_annotation.vcf \
    --output-file-format VCF \
    --data-sources-path ${launchDir}/funcotator_data/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38
    """
}