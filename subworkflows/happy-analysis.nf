include { DEEPVARIANT } from '../modules/deepvariant.nf'
include { HAPPY } from '../modules/happy.nf'


workflow ANALYSIS{
    take:
    ref
    fasta_index
    dict
    bench_bam
    bench_bam_idx
    wes_bed
    benchmark
    bench_idx
    bench_bed


    main:

    vcf_dv = Channel.empty()

    DEEPVARIANT(ref, bench_bam, bench_bam_idx, fasta_index, dict, wes_bed)
    vcf_dv = DEEPVARIANT.out.dv

    HAPPY(params.ref, params.fasta_index, params.benchmark, params.bench_idx, params.bench_bed, vcf_dv, params.wes_bed)
   
}