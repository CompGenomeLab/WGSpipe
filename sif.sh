#!/bin/bash
mkdir -p Singularity

# Build Singularity images from Docker containers
singularity build Singularity/bwa.sif docker://pegi3s/bwa
singularity build Singularity/deepvariant.sif docker://google/deepvariant:1.6.0
singularity build Singularity/ensembl.sif docker://ensemblorg/ensembl-vep
singularity build Singularity/fastp.sif docker://staphb/fastp
singularity build Singularity/fastqc.sif docker://staphb/fastqc
singularity build Singularity/freebayes.sif docker://staphb/freebayes
singularity build Singularity/gatk4.sif docker://broadinstitute/gatk
singularity build Singularity/happy.sif docker://jmcdani20/hap.py:v0.3.12
singularity build Singularity/multiqc.sif docker://staphb/multiqc
singularity build Singularity/sam.sif docker://dukegcb/bwa-samtools
singularity build Singularity/snpeff.sif docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_2
singularity build Singularity/picard.sif docker://broadinstitute/picard

#chmod +x sif.sh
#run ./sif.sh

