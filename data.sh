#!/bin/bash
#https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi;tab=live_object
#https://console.cloud.google.com/storage/browser/_details/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi;tab=live_object

#Make it executable chmod +x data.sh

#to run the script ./data.sh

# Define the data directory
data_dir="data"

# Download the Homo_sapiens_assembly38.fasta file using wget
wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz

wget -P "$data_dir" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

wget -P "$data_dir" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz

wget -P "$data_dir" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi

wget -P "$data_dir" ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

wget -P "$data_dir" ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

wget -P "$data_dir" "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list"

# Define the directory name
VEP_DATA_DIR="vep-data"

# Create the directory
mkdir -p $VEP_DATA_DIR

# Change to the directory
cd $VEP_DATA_DIR

# Download the VEP cache file
curl -O ftp://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz

# Extract the tar.gz file
tar xzf homo_sapiens_vep_110_GRCh38.tar.gz

# Optionally, remove the tar.gz file after extraction
rm homo_sapiens_vep_110_GRCh38.tar.gz

echo "VEP cache downloaded and extracted to $VEP_DATA_DIR"


