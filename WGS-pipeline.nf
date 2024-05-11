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
trimmomatic_adapters = params.genomes[params.genome]?.trimmomatic_adapters
trimmomatic_adapters_param = params.genomes[params.genome]?.trimmomatic_adapters_param
trimmomatic_window_len = params.genomes[params.genome]?.trimmomatic_window_len
trimmomatic_window_val = params.genomes[params.genome]?.trimmomatic_window_val
trimmomatic_min_len = params.genomes[params.genome]?.trimmomatic_min_len


params.benchmark = "${launchDir}/nf-core/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.bench_idx = "${launchDir}/nf-core/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
params.bench_bed = "${launchDir}/nf-core/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

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

    FASTQC(read_pairs_ch) // checks quality of the pairend reads fastq
    MULTIQC(FASTQC.out.fastqc_files) // checks quality of the pairend reads fastq
    
    BWA_INDEX(params.ref) //index reference for BWA
    TRIM_FASTP(read_pairs_ch) //trim paired end reads on default mode
    TRIMMOMATIC(read_pairs_ch, trimmomatic_adapters, trimmomatic_adapters_param, trimmomatic_window_len, trimmomatic_window_val, trimmomatic_min_len)
    BWA_MEM(params.ref, BWA_INDEX.out.index, TRIM_FASTP.out.trimmed) //
    SAM_CONVERTER(BWA_MEM.out.aligned) //converts to sam and sort the bam

    MARK_DEDUP(SAM_CONVERTER.out.bam, SAM_CONVERTER.out.bam_bai) //sorted bam file removed duplicates
    BASE_RECAP(MARK_DEDUP.out.dedup_bam, MARK_DEDUP.out.bai, params.ref, params.fasta_index, params.dict, params.site, params.site_idx, params.interval, params.dbsnp, params.dbsnp_idx, params.indel, params.indel_idx)
    APPLY_BQSR(MARK_DEDUP.out.dedup_bam, params.ref, BASE_RECAP.out.table, params.fasta_index, params.dict, params.interval)

    VARIANT_CALLING(params.ref, params.fasta_index, params.dict, APPLY_BQSR.out.applyed_bqsr_bam, APPLY_BQSR.out.bqsr_idx, params.dbsnp, params.dbsnp_idx, params.interval)
    vcf_htvc = VARIANT_CALLING.out.vcf_htvc
    vcf_dv = VARIANT_CALLING.out.vcf_dv
    vcf_fb = VARIANT_CALLING.out.vcf_fb
    vcf_all = VARIANT_CALLING.out.vcf_all

    HAPPY(params.ref, params.fasta_index, params.benchmark, params.bench_idx, params.bench_bed, vcf_dv)
    
    VAR_RECAL(params.ref, params.fasta_index, params.dict, vcf_htvc, params.dbsnp, params.thousandG, params.dbsnp_idx, params.thousandG_idx, params.hapmap, params.hapmap_idx, params.omni, params.omni_idx)
    APPLY_VQSR(params.ref,params.fasta_index, params.dict, vcf_htvc, VAR_RECAL.out.var_recal, VAR_RECAL.out.tranches, VAR_RECAL.out.var_recal_idx)
    VARIANT_FILTER(params.ref, params.fasta_index, params.dict, APPLY_VQSR.out.htvc_recalibrated, APPLY_VQSR.out.htvc_index_recalibrated)
    
    vcf_htvc = VARIANT_FILTER.out.htvc_filtered
  
    vcf_to_idx = Channel.empty()
    vcf_to_idx = vcf_to_idx.mix(vcf_htvc)
    vcf_to_idx = vcf_to_idx.mix(vcf_dv)
    vcf_to_idx = vcf_to_idx.mix(vcf_fb)

    GATK_INDEX(vcf_to_idx)
    vcf_index = GATK_INDEX.out.vcf_idx
    
    vcf_annotation = Channel.empty()
    vcf_annotation = vcf_annotation.mix(VARIANT_FILTER.out.htvc_filtered)
    vcf_annotation = vcf_annotation.mix(VARIANT_CALLING.out.vcf_dv)
    vcf_annotation = vcf_annotation.mix(VARIANT_CALLING.out.vcf_fb)

    ANNOTATION(params.ref, params.fasta_index, params.dict, vcf_annotation, vcf_index)
    
}

/*

nextflow run nf-core/sarek -profile singularity --input samplesheet.csv --outdir ./X-sarek --tools deepvariant,freebayes,haplotypecaller,snpeff,strelka,tiddit,manta,cnvkit,merge -r 3.4.0 -resume

include { SELECT_SNP_HTVC; SELECT_INDEL_HTVC } from './modules/select-variants-htvc.nf'
include { SELECT_SNP_DV; SELECT_INDEL_DV } from './modules/select-variants-dv.nf'
include { MERGE_VCFS_HTVC } from './modules/merge-vcf-htvc.nf'
include { MERGE_VCFS_DV } from './modules/merge-vcf-dv.nf'
include { FUNCOTATOR_ANNOTATION_HTVC } from './modules/funcotator-htvc.nf'
include { FUNCOTATOR_ANNOTATION_DV } from './modules/funcotator-dv.nf'
include { ANNOVAR_HTVC } from './modules/annovar-htvc.nf'
include { ANNOVAR_DV } from './modules/annovar-dv.nf'
include { SNPEFF_HTVC; SNPEFF_DV } from './modules/snpeff.nf'
    SELECT_SNP_HTVC(params.ref, params.fasta_index, params.dict, VARIANT_FILTER.out.htvc_filtered)
    SELECT_INDEL_HTVC(params.ref, params.fasta_index, params.dict, VARIANT_FILTER.out.htvc_filtered)

    SELECT_SNP_DV(params.ref, params.fasta_index, params.dict, vcf_dv)
    SELECT_INDEL_DV(params.ref, params.fasta_index, params.dict, vcf_dv)

    MERGE_VCFS_HTVC(SELECT_SNP_HTVC.out.htvc_snp, SELECT_INDEL_HTVC.out.htvc_indel)
    MERGE_VCFS_DV(SELECT_SNP_DV.out.dv_snp, SELECT_INDEL_DV.out.dv_indel)

    FUNCOTATOR_ANNOTATION_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    FUNCOTATOR_ANNOTATION_DV(params.ref, params.fasta_index, params.dict, MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)

    ANNOVAR_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    ANNOVAR_DV(params.ref, params.fasta_index, params.dict,  MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)

    SNPEFF_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx, snpeff_db)
    SNPEFF_DV(params.ref, params.fasta_index, params.dict, MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx, snpeff_db)
*/

/*
Reference genome options
  ascat_genome           : hg38
  ascat_alleles          : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_alleles_hg38.zip
  ascat_loci             : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_loci_hg38.zip
  ascat_loci_gc          : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/GC_G1000_hg38.zip
  ascat_loci_rt          : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/RT_G1000_hg38.zip
  bwa                    : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/
  bwamem2                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index/
  chr_dir                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes
  dbsnp                  : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz
  dbsnp_tbi              : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz.tbi
  dbsnp_vqsr             : --resource:dbsnp,known=false,training=true,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz
  dict                   : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict
  dragmap                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/dragmap/
  fasta                  : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
  fasta_fai              : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai
  germline_resource      : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz
  germline_resource_tbi  : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz.tbi
  known_indels           : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz
  known_indels_tbi       : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz.tbi
  known_indels_vqsr      : --resource:gatk,known=false,training=true,truth=true,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
  known_snps             : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz
  known_snps_tbi         : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz.tbi
  known_snps_vqsr        : --resource:1000G,known=false,training=true,truth=true,prior=10.0 1000G_omni2.5.hg38.vcf.gz
  mappability            : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/Control-FREEC/out100m2_hg38.gem
  ngscheckmate_bed       : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/NGSCheckMate/SNP_GRCh38_hg38_wChr.bed
  pon                    : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz
  pon_tbi                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi
  sentieon_dnascope_model: s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/Sentieon/SentieonDNAscopeModel1.1.model
  snpeff_db              : 105
  snpeff_genome          : GRCh38
  vep_genome             : GRCh38
  vep_species            : homo_sapiens
  vep_cache_version      : 110

nextflow run nf-core/sarek -r 3.4.0 -name test_hg002_3 -profile singularity -resume -params-file sarek/nf-params.json


snpEff \\
        -Xmx${avail_mem}M \\
        download ${genome}.${cache_version} \\
        -dataDir \${PWD}/snpeff_cache \\
        ${args}


    FUNCOTATOR_ANNOTATION_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    FUNCOTATOR_ANNOTATION_DV(params.ref, params.fasta_index, params.dict, MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)

    ANNOVAR_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    ANNOVAR_DV(params.ref, params.fasta_index, params.dict,  MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)

    SNPEFF_DB()
    snpeff_db = Channel.fromPath("\${PWD}/snpeff-db/*")

    SNPEFF_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx, snpeff_db)
    SNPEFF_DV(params.ref, params.fasta_index, params.dict, MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx, snpeff_db)
*/
