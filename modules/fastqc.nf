process FASTQC {
  tag "Fastqc process underway"

  publishDir "${params.outdir}/fastqc_reports/", mode: 'copy'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("*_fastqc.{zip,html}"), emit:fastqc_files

  script:
  """
  fastqc -t ${task.cpus} ${reads}
  """
}