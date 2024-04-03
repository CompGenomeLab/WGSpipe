process MERGE_VCFS_HTVC {

    tag "Merge HTVC Variant Files"

    publishDir "${params.outdir}/htvc-merge-variants", mode: 'copy'

    input:
    path snp
    path indel

    output:
    path "htvc_merged_variants.vcf", emit:htvc_merged
    path "*.vcf.idx", emit:htvc_merged_idx

    script:
    """
    gatk MergeVcfs -I ${snp} -I ${indel} -O htvc_merged_variants.vcf
    gatk IndexFeatureFile -I htvc_merged_variants.vcf
    """
}