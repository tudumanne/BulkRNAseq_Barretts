## Bulk RNA-seq pipeline for pathway analysis and biomarker discovery

This repository contains the data analysis pipeline used for the following project:

#### Transcriptomic Changes and Biomarkers in Barrett’s Metaplasia, Dysplasia and Cancer

[Add citation]

Data generation: 
The libraries were prepared using Illumina Stranded Total RNA Prep with Ribo-Zero Plus kit and sequenced on NovaSeq X with 2x150bp setting.

### Overview of the analysis workflow

Initial processing of sequencing data, including fastq read alignment and count estimations were run on a high-performance computing (HPC) system based on Linux. 

Downstream analyses were performed in R using RStudio on macOS Ventura 13.0.

The `hpc/` and `R/` folders contain additional information on the methods and the template scripts used to perform the analysis.

### Repository structure 
- `hpc/` - shell scripts for sequencing data processing 
- `R/` - R scripts for downstream analyses
- `annotation_files/` - custom annotation files
- `metadata/` - sample metadata

### Data avaiability 
Raw FASTQ files generated in this study have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEBxxxx (https://www.ebi.ac.uk/ena/browser/view/PRJEBxxxx).

Sample group information (`metadata/`) and custom annotation files (`annotation_files/`) are currently accessible, while processed data and raw data on ENA remain restricted until formal publication.

