process MULTIQC {
  tag "MULTIQC process underway"

  publishDir "${params.outdir}/multiqc/", mode: 'copy'
  
  input:
  tuple val(sample), path(reports)

  output:
  path("*.html")

  script:
  """
  multiqc ${reports}
  """
}