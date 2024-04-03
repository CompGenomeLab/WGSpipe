process MERGE_VCFS_DV {

    tag "Merge DV Variant Files"

    publishDir "${params.outdir}/dv-merge-variants", mode: 'copy'

    input:
    path snp
    path indel

    output:
    path "dv_merged_variants.vcf", emit:dv_merged
    path "*.vcf.idx", emit:dv_merged_idx

    script:
    """
    gatk MergeVcfs -I ${snp} -I ${indel} -O dv_merged_variants.vcf
    gatk IndexFeatureFile -I dv_merged_variants.vcf
    """
}