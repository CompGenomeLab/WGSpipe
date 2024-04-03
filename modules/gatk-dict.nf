process GATK_CREATE_DICTIONARY{
    tag "Creating dictionary of the reference"

    publishDir "${params.outdir}/gatk-ref-dict/", mode: 'copy'

    input:
    path ref

    output:
    path ("*.dict"), emit:ref_dict

    script:
    """
    gatk CreateSequenceDictionary -R ${ref}
    """
}
