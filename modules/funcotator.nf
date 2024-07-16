process FUNCOTATOR {

    tag "Annotating variants with Funcotator"

    publishDir "${params.outdir}/annotation/${merge_variant.baseName}_funcotator", mode: 'copy'

    input:
    path ref
    path fasta_index
    path dict
    path merge_variant
    path merge_variant_idx

    output:
    path "*.vcf"

    script:
    """
    gatk Funcotator \
    -R ${ref} \
    -V ${merge_variant}\
    -O ${merge_variant.baseName}_annotation.vcf \
    --output-file-format VCF \
    --data-sources-path ${launchDir}/funcotator_data/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38
    """
}