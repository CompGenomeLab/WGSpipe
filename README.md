# **WGS PIPELINE WITH GATK BEST PRACTICES**

This report is for training purposes on what steps should have been taken for GATK best practices. For constructing of the WGS pipeline GATK’s best practice workflow take into consideration with additional genomic analysis and annotation tools. The practices are written in Nextflow DSL2 workflow management language.

A whole genome sequence data should be picked from genomic databases such as 1000 Genomes Project or SRA database. Download data:

```jsx
Wget [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_1.fastq.gz)
Wget [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_2.fastq.g](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_2.fastq.gz)
```

Now we have our fastq files which corresponds to pair end reads of a WGS data.

Check the quality before getting to data preprocessing steps of the pipeline.

```jsx

  fastqc -t ${task.cpus} ${pairend_reads}
```

To evaluate your analysis to pass to further process go to this website below:

[Index of /projects/fastqc/Help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

# Step 1: Sorting and Alignment

Download a reference hg38 WGS data and other additional datas from google cloud provided by GATK:

> [https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0;tab=objects?pageState=("StorageObjectListTable":("f":"%255B%255D")](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)))&prefix=&forceOnObjectsSortingFiltering=false
> 

Index the reference genome using BWA-index tool:

```jsx
bwa index ${ref}"
```

- File types from bwa-index contains: `.amb, .ann, .bwt, .pac, .sa`

After indexing reference genome do the alignment using BWA-mem:

```jsx
bwa mem -t ${task.cpus} -R "@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}" ${ref} ${reads} > ${sample}_paired.sam
```

- @RG is an important command for further processes of the pipeline. It specifies the header of the read group which contains the identification, platform and sample name of the pair end read groups. Useful for interpertation

Converting and sorting the sam file that is generated through BWA-mem using sam tools.

```groovy
samtools view -bS ${sam_file} -o ${sam_file.baseName}_aligned.bam //convert
samtools sort -@ ${task.cpus} ${sam_file.baseName}_aligned.bam -o ${sam_file.baseName}_sorted_aligned.bam //sort
```

# **Step 2: Preprocessing of the Alignment file using GATK tools**

As following the GATK best practices these steps are important and necessary. These steps will allows us to improve the pipeline downstream analysis and generate accurate vcf files.

1. Marking duplicates and removal
2. Base recalibration
3. Haplotypecaller ( or a functionally same tool for variant call)

**MarkDuplicatesSpark**

Why is it neccessary to remove duplicates from the BAM files? — Removal of duplicates will help us to prevent PCR duplicate errors and allow us to proceed accurate variant calling in the end. 

```jsx
gatk MarkDuplicatesSpark \
-I ${bam_file} \ //aligned bam file
-O ${bam_file.baseName}_dedup.bam \ //output file name 
--remove-sequencing-duplicates \ //will remove duplicates from the file
--conf 'spark.executor.cores=${task.cpus}' // specify threads
```

**BaseRecalibratorSpark - ApplyBQSR**

Why is base recalibration important? — Base quality scores that have produced by sequencing machines can be prone to errors. Base recalibration and ApplyBQSR tool is able to correct these errors and improve variant calling.

You need to download a set of reference known sites data for Base Reclibration step. 

Known sites are needed for base recalibration to provide a reference dataset of variants that are believed to be true positives and for assessing the accuracy of base quality scores and identifying systematic errors in base calling.

```jsx
gatk BaseRecalibratorSpark \
-I ${dedup_bam} \
-R ${ref}  \
--known-sites ${sites} \
-O ${dedup_bam.baseName}_recal_data.table

gatk ApplyBQSR \
-I ${dedup_bam} \
-R ${ref}  \
--bqsr-recal-file ${table} \
-O ${dedup_bam.baseName}_recalibrated.bam

gatk --java-options "-Xmx4g" BuildBamIndex 
-I ${dedup_bam.baseName}_recalibrated.bam \
-O ${dedup_bam.baseName}_recalibrated.bam.bai //indexing the recalibrated file
```