# STR detection and reporting pipeline

## Introduction
This pipeline is designed to detect and report STRs using [ExpansionHunter](https://github.com/Illumina/ExpansionHunter), a powerful bioinformatics tool designed to detect repeat expansions from sequencing data.

## Pipeline Overview
The pipeline is structured to convert FASTQ files to BAM format, apply Expansion Hunter to detect specified loci, and finally execute a Python script for report generation.

## Pipeline Components
1. **FASTQ to BAM Conversion:**
   - The pipeline takes input in FASTQ format
   - It generates a SAM file
    - The SAM file is then converted to BAM format and then the BAM file is sorted

2. **Expansion Hunter Analysis:**
   - After the conversion, Expansion Hunter is utilized for the detection of specified loci
   - The loci are defined in the Expansion Hunter variant catalog

3. **Report Generation:**
   - Following the Expansion Hunter analysis, a Python script is executed to generate a comprehensive report

## Pipeline Parameters
1. **Input:**
    - FASTQ reads
    - FASTA reference file 
    - Expansion Hunter variant catalog

2. **Output:**
    - SAM file
    - BAM file
    - JSON file (Expansion Hunter output)
    - Report in PDF

## Running the Pipeline
```
nextflow run main.nf -c nextflow.config
```

### ** Side note:
Due to the size of the input files, all fastq file pairs and genome files are not included in this repository. 

Please upload the files in the *input* directory and replace the paths in the `nextflow.config` file with the appropriate paths to the input files.