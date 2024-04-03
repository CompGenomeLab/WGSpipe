process ANNOVAR_HTVC{
    tag "Annotating htvc variants with ANNOVAR"
    publishDir "${params.outdir}/annovar_htvc", mode: 'copy'

    input:
    path ref
    path fasta_index // Index of ref
    path dict // Dictionary of ref
    path htvc_variant // Input variant file (VCF format)
    path htvc_variant_idx

    output:
    path "*.hg38_multianno.vcf"
    path "*.hg38_multianno.txt"
    path "*.avinput"

    script:
    """
    table_annovar.pl ${htvc_variant} ${launchDir}/humandb/ \
        --buildver hg38 \
        --out ${htvc_variant.baseName}_annovar \
        --remove \
        --protocol refGene,gnomad_genome,gnomad_exome,avsnp150,clinvar_20221231 \
        --operation g,f,f,f,f \
        --nastring . \
        --vcfinput \
        --thread ${task.cpus}
    """
}