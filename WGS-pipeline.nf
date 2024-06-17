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
params.trimmomatic_adapters = params.genomes[params.genome]?.trimmomatic_adapters
params.trimmomatic_adapters_param = params.genomes[params.genome]?.trimmomatic_adapters_param
params.trimmomatic_window_len = params.genomes[params.genome]?.trimmomatic_window_len
params.trimmomatic_window_val = params.genomes[params.genome]?.trimmomatic_window_val
params.trimmomatic_min_len = params.genomes[params.genome]?.trimmomatic_min_len
params.tools = params.genomes[params.genome]?.tools
params.step = params.genomes[params.genome]?.step
params.bam_file = params.bam_file[params.genome]?.bam_file
params.bam_file_idx = params.bam_file_idx[params.genome]?.bam_file_idx
params.vcf_idx = params.vcf_idx[params.genome]?.vcf_idx
params.vcf_annotation = params.vcf_annotation[params.genome]?.vcf_annotation


params.bam = "${launchDir}/nf-core/input/HG003.novaseq.wes_idt.100x.dedup.bam"
params.bam_bai = "${launchDir}/nf-core/input/HG003.novaseq.wes_idt.100x.dedup.bam"
params.benchmark = "${launchDir}/nf-core/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.bench_idx = "${launchDir}/nf-core/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
params.bench_bed = "${launchDir}/nf-core/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
params.wes_bed = "${launchDir}/nf-core/benchmark/idt_capture.grch38.bed"

read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true) //paired end check

def helpMessage() {
    log.info"""
    =========================================
     BaharSevgin/Multiple-Variant-Call-Pipeline
    =========================================
    Usage:
    
    The typical command for running the pipeline is as follows:
    
    nextflow run WGS-pipeline.nf --reads '*_R{1,2}.fastq.gz' 
    Mandatory arguments:
      --reads                      Path to input data (must be surrounded with quotes).
      --ref                        Path to human Fasta reference
    References
      --fasta_index                If reference index files exist.
      --dict                       If reference dictionary files exist.
    Options:
    --step                         
    
    Trimming options
      --trimmomatic_adapters       Adapters index for adapter removal
      --trimmomatic_adapters_param Trimming parameters for adapters. <seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. Default 2:30:10
      --trimmomatic_window_len     Window size. Default 4
      --trimmomatic_window_val     Window average quality requiered. Default 20
      --trimmomatic_min_len        Minimum length of reads
    Other options:
      --outdir                     The output directory where the results will be saved
    """.stripIndent()
}

params.help = false
if (params.help) {
    helpMessage()
    exit 0
}

log.info """\
      WGS PİPELİNE
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
include { MARK_DEDUP } from './modules/mark-dedup.nf'
include { BASE_RECAP } from './modules/base-recalibrate.nf'
include { APPLY_BQSR } from './modules/apply-bqsr.nf'

include {VARIANT_CALLING} from './subworkflows/variant-call-all-germline.nf'

include { HAPPY } from './modules/happy.nf'
include { VAR_RECAL } from './modules/variant-recalibrate.nf'
include { APPLY_VQSR } from './modules/apply-vqsr.nf'
include { VARIANT_FILTER } from './modules/variant-filteration.nf'
include { GATK_INDEX } from './modules/gatk4-index.nf'

include {ANNOTATION} from './subworkflows/annotation-all.nf'
include {ANALYSIS} from './subworkflows/happy-analysis.nf'


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


    if(params.step.split(,).contains('preprocessing')){
        FASTQC(read_pairs_ch) // checks quality of the pairend reads fastq
        MULTIQC(FASTQC.out.fastqc_files) // checks quality of the pairend reads fastq
        BWA_INDEX(params.ref) //index reference for BWA
        
        TRIM_FASTP(read_pairs_ch) //trim paired end reads on default mode
        //TRIMMOMATIC(read_pairs_ch, params.trimmomatic_adapters, params.trimmomatic_adapters_param, params.trimmomatic_window_len, params.trimmomatic_window_val, params.trimmomatic_min_len)
        
        BWA_MEM(params.ref, BWA_INDEX.out.index, TRIM_FASTP.out.trimmed) 
        SAM_CONVERTER(BWA_MEM.out.aligned) //converts to sam and sort the bam

        MARK_DEDUP(SAM_CONVERTER.out.bam, SAM_CONVERTER.out.bam_bai) //sorted bam file removed duplicates
        BASE_RECAP(MARK_DEDUP.out.dedup_bam, MARK_DEDUP.out.bai, params.ref, params.fasta_index, params.dict, params.site, params.site_idx, params.interval, params.dbsnp, params.dbsnp_idx, params.indel, params.indel_idx)
        APPLY_BQSR(MARK_DEDUP.out.dedup_bam, params.ref, BASE_RECAP.out.table, params.fasta_index, params.dict, params.interval)
        
        params.bam_file = APPLY_BQSR.out.applyed_bqsr_bam
        params.bam_file_idx = APPLY_BQSR.out.bqsr_idx
    
    }
    
    vcf_index_ch = Channel.empty()
    vcf_annotation_ch = Channel.empty()
    
    if(params.step.split(,).contains('variant_calling')){
      if(!params.bam_file || !params.bam_file_idx){ //ref file check
        exit 1, "BAM file of a genome not specified! Please, provide --bam_file and --bam_file_idx"
      }
      if (params.tools.split(',').contains('haplotypecaller') || params.tools.split(',').contains('deepvariant') || params.tools.split(',').contains('freebayes')){
          VARIANT_CALLING(params.tools, params.ref, params.fasta_index, params.dict, bam_file, bam_file_idx, params.dbsnp, params.dbsnp_idx, params.interval)
          vcf_htvc = VARIANT_CALLING.out.vcf_htvc
          vcf_dv = VARIANT_CALLING.out.vcf_dv
          vcf_fb = VARIANT_CALLING.out.vcf_fb
          vcf_all = VARIANT_CALLING.out.vcf_all

          HAPPY(params.ref, params.fasta_index, params.benchmark, params.bench_idx, params.bench_bed, vcf_dv, params.wes_bed)
          
          VAR_RECAL(params.ref, params.fasta_index, params.dict, vcf_htvc, params.dbsnp, params.thousandG, params.dbsnp_idx, params.thousandG_idx, params.hapmap, params.hapmap_idx, params.omni, params.omni_idx)
          APPLY_VQSR(params.ref,params.fasta_index, params.dict, vcf_htvc, VAR_RECAL.out.var_recal, VAR_RECAL.out.tranches, VAR_RECAL.out.var_recal_idx)
          VARIANT_FILTER(params.ref, params.fasta_index, params.dict, APPLY_VQSR.out.htvc_recalibrated, APPLY_VQSR.out.htvc_index_recalibrated)
          
          vcf_htvc = VARIANT_FILTER.out.htvc_filtered
  
          vcf_to_idx = vcf_to_idx.mix(vcf_htvc)
          vcf_to_idx = vcf_to_idx.mix(vcf_dv)
          vcf_to_idx = vcf_to_idx.mix(vcf_fb)

          GATK_INDEX(vcf_to_idx)
          vcf_index_ch = GATK_INDEX.out.vcf_idx
          
          vcf_annotation_ch = vcf_annotation.mix(VARIANT_FILTER.out.htvc_filtered)
          vcf_annotation_ch = vcf_annotation.mix(VARIANT_CALLING.out.vcf_dv)
          vcf_annotation_ch = vcf_annotation.mix(VARIANT_CALLING.out.vcf_fb)
      }

    }

    
    if (!params.steps.contains('variant_calling') && !params.steps.contains('preprocessing')) {  
      vcf_annotation_ch = Channel.fromFilePairs(params.vcf_annotation, checkIfExists:true)
      vcf_index_ch = Channel.fromFilePairs(params.vcf_idx, checkIfExists:true)
    }
    
    if(params.step.split(,).contains('annotate')){
      if(!params.vcf_idx || !params.vcf_annotation){
         exit 1, "VCF file of a genome not specified for annotation! Please, provide --vcf_annotation and --vcf_idx"
      }
          

      if (params.tools.split(',').contains('funcotator') || params.tools.split(',').contains('annovar') || params.tools.split(',').contains('snpeff')){
         ANNOTATION(params.tools, params.ref, params.fasta_index, params.dict, vcf_annotation, vcf_index, params.snpeff_db)
      }
    }
  
}


/*
NOTES

salloc -A users -p long_mdbf --qos=long_mdbf $SHELL
scontrol show job 2075|grep NodeList
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/ 
https://console.cloud.google.com/storage/browser/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))

nextflow run nf-core/sarek -profile singularity --input samplesheet.csv --outdir ./X-sarek --tools deepvariant,freebayes,haplotypecaller,snpeff,strelka,tiddit,manta,cnvkit,merge -r 3.4.0 -resume
*/
