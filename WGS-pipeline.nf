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

params.benchmark = "${launchDir}/new_data/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.bench_idx = "${launchDir}/new_data/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
params.bench_bed = "${launchDir}/new_data/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

log.info """\
      WGS PİPELİNE
================================
Reference        : ${params.ref}
Reads            : ${params.reads}
Output-folder    : ${params.outdir}/
Known-sites      : ${params.site}
"""

if (params.ref) { 

  dbsnp        = file( params.dbsnp )
  site         = file( params.site )
  indel        = file( params.indel )
  hapmap       = file( params.hapmap )
  thousandG    = file( params.thousandG )
  omni         = file( params.omni )
  dbsnp_idx    = file( params.dbsnp_index )
  site_idx     = file( params.site_idx )
  indel_idx    = file( params.indel_idx )
  hapmap_idx   = file( params.hapmap_idx )
  thousandG_idx  = file( params.phase_idx )
  omni_idx     = file( params.omni_idx )

 }  else { exit 1, 'ERROR: Goole cloud reference directory for GATK variants calling is not specified' }


include { SAM_INDEX_REF_FASTA } from './modules/sam-index.nf'
include { GATK_CREATE_DICTIONARY } from './modules/gatk-dict.nf'
include { FASTQC } from './modules/fastqc.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { BWA_INDEX } from './modules/bwa-index.nf'
include { TRIM_FASTP } from './modules/fastp.nf'
include { BWA_MEM } from './modules/bwa-mem.nf'
include { SAM_CONVERTER } from './modules/samtools.nf'
include { MARK_DEDUP } from './modules/mark-dedup.nf'
include { BASE_RECAP } from './modules/base-recalibrate.nf'
include { APPLY_BQSR } from './modules/apply-bqsr.nf'
include { HAPLOTYPECALLER } from './modules/haplotypecaller.nf'
include { DEEPVARIANT } from './modules/deepvariant.nf'
include { HAPPY } from './modules/happy.nf'
include { VAR_RECAL } from './modules/variant-recalibrate.nf'
include { APPLY_VQSR } from './modules/apply-vqsr.nf'
include { VARIANT_FILTER } from './modules/variant-filteration.nf'
include { SELECT_SNP_HTVC, SELECT_INDEL_HTVC } from './modules/select-variants-htvc.nf'
include { SELECT_SNP_DV, SELECT_INDEL_DV } from './modules/select-variants-dv.nf'
include { MERGE_VCFS_HTVC } from './modules/merge-vcf-htvc.nf'
include { MERGE_VCFS_DV } from './modules/merge-vcf-dv.nf'
include { FUNCOTATOR_ANNOTATION_HTVC } from './modules/funcotator-htvc.nf'
include { FUNCOTATOR_ANNOTATION_DV } from './modules/funcotator-dv.nf'
include { ANNOVAR_HTVC } from './modules/annovar-htvc.nf'
include { ANNOVAR_DV } from './modules/annovar-dv.nf'

workflow {

    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)
    read_pairs_ch.view()

    if(params.fasta_index!){
        index = SAM_INDEX_REF_FASTA(params.ref)
        params.fasta_index = index
    }

    if(params.dict!){
       dict = GATK_CREATE_DICTIONARY(params.ref)
       params.dict = dict
    }

    FASTQC(read_pairs_ch)
    FASTQC.out.fastqc_files.view()
    MULTIQC(FASTQC.out.fastqc_files)
    
    BWA_INDEX(params.ref)
    TRIM_FASTP(read_pairs_ch)
    BWA_MEM(params.ref, BWA_INDEX.out.index, TRIM_FASTP.out.trimmed)    
    SAM_CONVERTER(BWA_MEM.out.aligned)

    MARK_DEDUP(SAM_CONVERTER.out.bam, SAM_CONVERTER.out.bam_bai)
    BASE_RECAP(MARK_DEDUP.out.dedup_bam, MARK_DEDUP.out.bai, params.ref, params.fasta_index, params.dict, params.site, params.site_idx, params.interval, params.dbsnp, params.dbsnp_idx, params.indel, params.indel_idx)
    APPLY_BQSR(MARK_DEDUP.out.dedup_bam, params.ref, BASE_RECAP.out.table, params.fasta_index, params.dict, params.interval)

    HAPLOTYPECALLER(params.ref, APPLY_BQSR.out.applyed_bqsr_bam, APPLY_BQSR.out.bqsr_idx, params.fasta_index, params.dict,params.dbsnp, params.dbsnp_idx, params.interval)
    DEEPVARIANT(params.ref, APPLY_BQSR.out.applyed_bqsr_bam, APPLY_BQSR.out.bqsr_idx, params.fasta_index, params.dict)

    HAPPY(params.ref, params.fasta_index, params.benchmark, params.bench_idx, params.bench_bed, DEEPVARIANT.out.dv)
    
    VAR_RECAL(params.ref, params.fasta_index, params.dict, HAPLOTYPECALLER.out.htvc, params.dbsnp, params.thousandG, params.dbsnp_idx, params.thousandG_idx, params.hapmap, params.hapmap_idx, params.omni, params.omni_idx)
    APPLY_VQSR(params.ref,params.fasta_index, params.dict, HAPLOTYPECALLER.out.htvc, VAR_RECAL.out.var_recal, VAR_RECAL.out.tranches, VAR_RECAL.out.var_recal_idx)
    VARIANT_FILTER(params.ref, params.fasta_index, params.dict, APPLY_VQSR.out.htvc_recalibrated, APPLY_VQSR.out.htvc_index_recalibrated)

    SELECT_SNP_HTVC(params.ref, params.fasta_index, params.dict, VARIANT_FILTER.out.htvc_filtered)
    SELECT_INDEL_HTVC(params.ref, params.fasta_index, params.dict, VARIANT_FILTER.out.htvc_filtered)

    SELECT_SNP_DV(params.ref, params.fasta_index, params.dict, DEEPVARIANT.out.dv)
    SELECT_INDEL_DV(params.ref, params.fasta_index, params.dict, DEEPVARIANT.out.dv)

    MERGE_VCFS_HTVC(SELECT_SNP_HTVC.out.htvc_snp, SELECT_INDEL_HTVC.out.htvc_indel)
    MERGE_VCFS_DV(SELECT_SNP_DV.out.dv_snp, SELECT_INDEL_DV.out.dv_indel)

    FUNCOTATOR_ANNOTATION_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    FUNCOTATOR_ANNOTATION_DV(params.ref, params.fasta_index, params.dict, MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)

    ANNOVAR_HTVC(params.ref, params.fasta_index, params.dict, MERGE_VCFS_HTVC.out.htvc_merged, MERGE_VCFS_HTVC.out.htvc_merged_idx)
    ANNOVAR_DV(params.ref, params.fasta_index, params.dict,  MERGE_VCFS_DV.out.dv_merged, MERGE_VCFS_DV.out.dv_merged_idx)
}

