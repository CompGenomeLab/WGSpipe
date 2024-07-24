params.genome = 'hg38' 

params.ref = params.genomes[params.genome]?.ref
params.reads = params.genomes[params.genome]?.reads
params.site = params.genomes[params.genome]?.site
params.outdir = params.genomes[params.genome]?.outdir
params.dict = params.genomes[params.genome]?.dict
params.fasta_index = params.genomes[params.genome]?.fasta_index
params.site_idx = params.genomes[params.genome]?.site_idx
params.dbsnp = params.genomes[params.genome]?.dbsnp
params.thousandG = params.genomes[params.genome]?.thousandG
params.dbsnp_idx = params.genomes[params.genome]?.dbsnp_idx
params.thousandG_idx = params.genomes[params.genome]?.thousandG_idx
params.hapmap = params.genomes[params.genome]?.hapmap
params.hapmap_idx = params.genomes[params.genome]?.hapmap_idx
params.omni = params.genomes[params.genome]?.omni
params.omni_idx = params.genomes[params.genome]?.omni_idx
params.interval = params.genomes[params.genome]?.interval
params.indel = params.genomes[params.genome]?.indel
params.indel_idx = params.genomes[params.genome]?.indel_idx
params.snpeff_db = params.genomes[params.genome]?.snpeff_db 
params.snpeff_cache = params.genomes[params.genome]?.snpeff_cache  
params.tools = params.genomes[params.genome]?.tools
params.step = params.genomes[params.genome]?.step
params.bam_file = params.genomes[params.genome]?.bam_file
params.bam_file_idx = params.genomes[params.genome]?.bam_file_idx
params.vep_cache_version = params.genomes[params.genome]?.vep_cache_version
params.vep_genome = params.genomes[params.genome]?.vep_genome            
params.vep_species = params.genomes[params.genome]?.vep_species 

params.benchmark = "${launchDir}/nf-core/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.bench_idx = "${launchDir}/nf-core/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
params.bench_bed = "${launchDir}/nf-core/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"


read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true) //paired end check

def helpMessage() {
    log.info"""
    =========================================
     BaharSevgin/Multiple-Variant-Call-Pipeline
    =========================================
    Usage:
    
    The typical command for running the pipeline is as follows:
    
    nextflow run WGS-pipeline.nf --reads '*_R{1,2}.fastq.gz' --ref "*.fasta"
    Mandatory arguments:
      --reads                      Path to input data (must be surrounded with quotes).
      --ref                        Path to human Fasta reference
    References
      --fasta_index                If reference index files exist.
      --dict                       If reference dictionary files exists.
    Other options:
      --outdir                     The output directory where the results will be saved
      --step                       Define the steps you wanted to use. Default all included
      --tools                      Define the tools you wanted to use in the steps of variant calling and annotation. Default all included
    """.stripIndent()
}

params.help = false
if (params.help) {
    helpMessage()
    exit 0
}

log.info """\
 Multiple-Variant-Call-Pipeline
================================
Reference        : ${params.ref}
Reads            : ${params.reads}
Output-folder    : ${params.outdir}/
Known-sites      : ${params.site}
"""

include { SAM_INDEX_REF_FASTA } from './modules/sam-index.nf'
include { GATK_CREATE_DICTIONARY } from './modules/gatk-dict.nf'
include { FASTQC } from './modules/fastqc.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { BWA_INDEX } from './modules/bwa-index.nf'
include { TRIM_FASTP } from './modules/fastp.nf'
include { TRIMMOMATIC } from './modules/trimmomatic.nf'
include { BWA_MEM } from './modules/bwa-mem.nf'
include { SAM_CONVERTER } from './modules/samtools.nf'
include { MARK_DEDUP_PICARD } from './modules/mark-dedup-picard.nf'
include { BASE_RECAP_PICARD } from './modules/recalibration-picard.nf'
include { APPLY_BQSR_PICARD } from './modules/apply_bqsr_picard.nf'
include {VARIANT_CALLING} from './subworkflows/variant-call-all-germline.nf'
include { HAPPY } from './modules/happy.nf'
include { VAR_RECAL } from './modules/variant-recalibrate.nf'
include { APPLY_VQSR } from './modules/apply-vqsr.nf'
include { VARIANT_FILTER } from './modules/variant-filteration.nf'
include { GATK_INDEX } from './modules/gatk4-index.nf'
include { SNPEFF } from './modules/snpeff.nf'
include { ENSEMBL_VEP } from './modules/enesembl.nf'

workflow {
    if(!params.ref){ //ref file check
      exit 1, "Reference genome not specified! Please, provide --ref"
    }
    if(!params.reads){ //ref file check
      exit 1, "Pairedend reads of a genome not specified! Please, provide --reads"
    }

    if(params.ref != "${launchDir}/data/Homo_sapiens_assembly38.fasta"){ //if not default generate index and dict
        index = SAM_INDEX_REF_FASTA(params.ref)
        params.fasta_index = index

        dict = GATK_CREATE_DICTIONARY(params.ref)
        params.dict = dict
    }

    final_bam_ch = Channel.empty() //when preprocessing finishes
    final_bam_idx_ch = Channel.empty() //when preprocessing finishes

    if(params.step.split(',').contains('preprocessing')){
        FASTQC(read_pairs_ch) // checks quality of the pairend reads fastq
        MULTIQC(FASTQC.out.fastqc_files) // checks quality of the pairend reads fastq
        BWA_INDEX(params.ref) //index reference for BWA
        
        TRIM_FASTP(read_pairs_ch) //trim paired end reads on default mode
        
        BWA_MEM(params.ref, BWA_INDEX.out.index, TRIM_FASTP.out.trimmed) 
        SAM_CONVERTER(BWA_MEM.out.aligned) //converts to sam and sort the bam

        MARK_DEDUP_PICARD(SAM_CONVERTER.out.bam, SAM_CONVERTER.out.bam_bai) //sorted bam file removed duplicates
        BASE_RECAP_PICARD(MARK_DEDUP_PICARD.out.dedup_bam, MARK_DEDUP_PICARD.out.bai, params.ref, params.fasta_index, params.dict, params.site, params.site_idx, params.interval, params.dbsnp, params.dbsnp_idx, params.indel, params.indel_idx)
        APPLY_BQSR_PICARD(MARK_DEDUP_PICARD.out.dedup_bam, params.ref, BASE_RECAP_PICARD.out.table, params.fasta_index, params.dict, params.interval)
        
        final_bam_ch = APPLY_BQSR_PICARD.out.applyed_bqsr_bam
        final_bam_idx_ch = APPLY_BQSR_PICARD.out.bqsr_idx
    }
    else if (params.step.split(',').contains('variant_calling')) {
        // If preprocessing is skipped, use input BAM file
        if (!params.bam_file || !params.bam_file_idx) {
            exit 1, "BAM file of a genome not specified! Please, provide --bam_file and --bam_file_idx"
        }
        final_bam_ch = Channel.fromPath(params.bam_file)
        final_bam_idx_ch = Channel.fromPath(params.bam_file_idx)
    }

    vcf_index_ch = Channel.empty() //creating index after vcf 
    vcf_annotation_ch = Channel.empty() //creating a mixed channel for snpeff

    if(params.step.split(',').contains('variant_calling')){

        if (params.tools.split(',').contains('haplotypecaller') || params.tools.split(',').contains('deepvariant') || params.tools.split(',').contains('freebayes')){

          VARIANT_CALLING(params.tools, params.ref, params.fasta_index, params.dict, final_bam_ch, final_bam_idx_ch, params.dbsnp, params.dbsnp_idx, params.interval)
          vcf_htvc = VARIANT_CALLING.out.vcf_htvc
          vcf_dv = VARIANT_CALLING.out.vcf_dv
          vcf_fb = VARIANT_CALLING.out.vcf_fb

          HAPPY(params.ref, params.fasta_index, params.benchmark, params.bench_idx, params.bench_bed, vcf_dv)
          
          VAR_RECAL(params.ref, params.fasta_index, params.dict, vcf_htvc, params.dbsnp, params.thousandG, params.dbsnp_idx, params.thousandG_idx, params.hapmap, params.hapmap_idx, params.omni, params.omni_idx)
          APPLY_VQSR(params.ref,params.fasta_index, params.dict, vcf_htvc, VAR_RECAL.out.var_recal, VAR_RECAL.out.tranches, VAR_RECAL.out.var_recal_idx)
          VARIANT_FILTER(params.ref, params.fasta_index, params.dict, APPLY_VQSR.out.htvc_recalibrated, APPLY_VQSR.out.htvc_index_recalibrated)
          
          vcf_htvc = VARIANT_FILTER.out.htvc_filtered

          vcf_annotation_ch = vcf_annotation_ch.mix(VARIANT_FILTER.out.htvc_filtered)
          vcf_annotation_ch = vcf_annotation_ch.mix(VARIANT_CALLING.out.vcf_dv)
          vcf_annotation_ch = vcf_annotation_ch.mix(VARIANT_CALLING.out.vcf_fb)
  
          vcf_to_idx = Channel.empty() //vcf files that needs to be indexed 
          vcf_to_idx = vcf_to_idx.mix(vcf_htvc)
          vcf_to_idx = vcf_to_idx.mix(vcf_dv)
          vcf_to_idx = vcf_to_idx.mix(vcf_fb)

          GATK_INDEX(vcf_to_idx) //creating index of the vcf files
          vcf_index_ch = GATK_INDEX.out.vcf_idx //idx files transferred to a index_ch
              
        }
    }

    if(params.step.split(',').contains('annotate'))
    {
        SNPEFF(params.ref, params.fasta_index, params.dict, vcf_annotation_ch, vcf_index_ch, params.snpeff_db)
    }


}

