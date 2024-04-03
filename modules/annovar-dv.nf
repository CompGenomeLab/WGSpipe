process ANNOVAR_DV {
    tag "Annotating dv variants with ANNOVAR" 
    publishDir "${params.outdir}/annovar_dv", mode: 'copy'

    input:
    path ref
    path fasta_index // Index of ref
    path dict // Dictionary of ref
    path dv_variant // Input variant file (VCF format)
    path dv_variant_idx

    output:
    path "*.hg38_multianno.vcf"
    path "*.hg38_multianno.txt"
    path "*.avinput"

    script:
    """
    table_annovar.pl ${dv_variant} ${launchDir}/humandb/ \
        --buildver hg38 \
        --out ${dv_variant.baseName}_annovar \
        --remove \
        --protocol refGene,gnomad_genome,gnomad_exome,avsnp150,clinvar_20221231 \
        --operation g,f,f,f,f \
        --nastring . \
        --vcfinput \
        --thread ${task.cpus}
    """
}