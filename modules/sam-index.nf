process SAM_INDEX_REF_FASTA {
    tag "Creating index of the reference"

    publishDir "${params.outdir}/sam-ref-idx/", mode: 'copy'

	input:
	path ref

	output:
	path ("*.fai"),emit: genome_idx
	script:
	"""
	samtools faidx ${ref}
	"""
}