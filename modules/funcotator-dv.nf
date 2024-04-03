process FUNCOTATOR_ANNOTATION_DV {

    tag "Funcotator annotation on dv"

    publishDir "${params.outdir}/dv_funcotator_annotation", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx

    output:
    path "dv_funcotator_annotation.vcf"

    script:
    """
    gatk Funcotator \
    -R ${ref} \
    -V ${merge_variant}\
    -O dv_funcotator_annotation.vcf \
    --output-file-format VCF \
    --data-sources-path ${launchDir}/funcotator_data/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38
    """
}