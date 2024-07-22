#!/bin/bash

# Configure conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install Singularity
conda install conda-forge::singularity 

# Create a directory for Singularity images
mkdir -p Singularity

# Build Singularity images from Docker containers
singularity build Singularity/annovar.sif docker://bioinfochrustrasbourg/annovar
singularity build Singularity/bwa.sif docker://pegi3s/bwa
singularity build Singularity/deepvariant.sif docker://google/deepvariant:1.6.0
singularity build Singularity/fastp.sif docker://staphb/fastp
singularity build Singularity/fastqc.sif docker://staphb/fastqc
singularity build Singularity/freebayes.sif docker://staphb/freebayes
singularity build Singularity/gatk4.sif docker://broadinstitute/gatk
singularity build Singularity/happy.sif docker://jmcdani20/hap.py:v0.3.12
singularity build Singularity/multiqc.sif docker://staphb/multiqc
singularity build Singularity/sam.sif docker://dukegcb/bwa-samtools
singularity build Singularity/snpeff.sif docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_2

#chmod +x sif.sh
#run ./install_singularity_tools.sh

