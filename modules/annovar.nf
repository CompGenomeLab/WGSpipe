process ANNOVAR{
    tag "Annotating variants with ANNOVAR"
    publishDir "${params.outdir}/annovar_annotation", mode: 'copy'

    input:
    path ref
    path fasta_index // Index of ref
    path dict // Dictionary of ref
    path variant // Input variant file (VCF format)
    path variant_idx

    output:
    path "*.hg38_multianno.vcf"
    path "*.hg38_multianno.txt"
    path "*.avinput"

    script:
    """
    table_annovar.pl ${variant} ${launchDir}/humandb/ \
        --buildver hg38 \
        --out ${variant.baseName}_annovar \
        --remove \
        --protocol refGene,gnomad_genome,gnomad_exome,avsnp150,clinvar_20221231 \
        --operation g,f,f,f,f \
        --nastring . \
        --vcfinput \
        --thread ${task.cpus}
    """
}