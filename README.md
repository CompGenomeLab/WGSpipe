# **WGS PIPELINE WITH GATK BEST PRACTICES**

This report is for training purposes on what steps should have been taken for GATK best practices. For constructing of the WGS pipeline GATK’s best practice workflow take into consideration with additional genomic analysis and annotation tools. The practices are written in Nextflow DSL2 workflow management language.

A whole genome sequence data should be picked from genomic databases such as 1000 Genomes Project or SRA database. Download data:

```jsx
Wget [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_1.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_1.fastq.gz)
Wget [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_2.fastq.g](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/SRR794275/SRR794275_2.fastq.gz)
```

Now we have our fastq files which corresponds to pair end reads of a WGS data.

For the needed data resources the google bucket data set is used which provided by GATK:
https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

Check the quality before getting to data preprocessing steps of the pipeline.

```jsx

  fastqc -t ${task.cpus} ${pairend_reads}
```

To evaluate your analysis to pass to further process go to this website below:

[Index of /projects/fastqc/Help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

To check quality in a more detailed manner **multiqc** is used.

```jsx
multiqc ${directory or files}
```

Trimming step is proceed with **fastp** tool**.** Commonly default parameters were used, but can add specific parameters. Fastp have single end and paired end options.

```jsx
    fastp --in1 ${reads[0]} --in2 ${reads[1]} \
    --thread ${task.cpus} \
    --out1 ${sample}_trimmed_1.fastq --out2 ${sample}_trimmed_2.fastq \
    -j ${sample}_fastp.json -h ${sample}_fastp.html
```

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

**HaplotypeCaller — Variant Call tool**

Particularly made for germline variations found in WES or WGS data. Generates precisely accurate vcf files. Uses joint calling, handling many samples simultaneously.

```jsx
gatk HaplotypeCaller \
-R ${ref} \
-I ${applyed_bqsr_bam} \
-O htvc_variants.vcf
```

There are other tools that could be used in WGS or WES pipelines for variant calling. One of the major example is **DeepVariant.** Deepvariant ****aims to accurately identify variants in WGS and WES data and do its own variant filtration and recalibration.

```jsx
/opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${ref} \
        --reads=${applyed_bqsr_bam} \
        --output_vcf=dv_variants.vcf \
        --num_shards=${task.cpus}
```

# **Step 3: Preparing Variant File for Annotation**

There are 2 major steps for filtration and recalibration:

1. Variant Recalibrator & Apply VQSR
2. Variant Filtration 

For recalibration it needs resources to recalibrator with a truth data set which can be retrieved from GATK google cloud database. These truth datasets includes 1000GP, HapMap, dbSNP, OMNI projects and more.

```jsx
 gatk VariantRecalibrator \
            -R ${ref} \
            -V ${variants} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousandG} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode BOTH \
            -O ${variants.baseName}_output.recal \
            --tranches-file ${variants.baseName}_output.tranches

    gatk ApplyVQSR \
        -R ${ref} \
        -V ${variants}\
        --recal-file ${var_recal} \
        --tranches-file ${tranches} \
        -mode BOTH \
        -ts-filter-level 99.0 \
        -O ${variants.baseName}_recalibrated_variant.vcf
```

Variant Filtration going to be done after recalibrator process with proper parameters that is accepted by the GATK best practices.

```jsx
gatk VariantFiltration \
            -R ${ref} \
            -V ${recal_variants} \
            -O htvc_variant_filtered.vcf \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 60.0" \
            -filter-name "MQ_filter" -filter "MQ < 40.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
```

After processing the vcf file variants will be selected with Select Variants tool provided by GATK.

```jsx
gatk SelectVariants -R ${ref} -V ${variants} \
--select-type-to-include SNP \ //SNP, INDEL, or SNV
-O htvc_snp_filtered_variants.vcf
```

# STEP 4: ANNOTATING VARIANT FILE

Annotation can be done using various annotation tools: Annovar, Funcotator, or SnpEff. For this pipeline Annovar and Funcotator were included for comparison of accuracy of the annotation tools. 

For both tool various truth data sets were included. You can change them according to your purpose. 

**Funcotator:**

```jsx
  gatk Funcotator \
    -R ${reference} \
    -V ${variant}\
    -O annotation.vcf \
    --output-file-format VCF \
    --data-sources-path ${launchDir}/funcotator_data/funcotator_dataSources.v1.8.hg38.20230908g \
    --ref-version hg38
```

**ANNOVAR:**

```jsx
  table_annovar.pl ${variant} ${launchDir}/humandb/ \
        --buildver hg38 \
        --out annovar \
        --remove \
        --protocol refGene,gnomad_genome,gnomad_exome,avsnp150,clinvar_20221231 \
        --operation g,f,f,f,f \
        --nastring . \
        --vcfinput \
        --thread ${task.cpus}
```