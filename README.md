# RNA-seq Analysis Pipeline for Pathway Analysis and Biomarker Discovery in Barrett's Neoplasia

This repository contains the data analysis pipeline used for the analysis of bulk RNA-seq data described in:
[Add citation]

### Overview of the analysis workflow

Sequencing data alignment and read count estimations were run on a high-performance computing (HPC) system based on Linux. 
R based analyses were carried out using RStudio run on MacOS Ventura 13.0.

The 'hpc' and 'R' folders contain more information on how the analyses were run and template scripts used.

### Workflow

``` mermaid
graph TD; 
    A[Raw FASTQ files]-->B[FastQC quality check];
    A-->C[STAR alignment];
    A-->D[Salmon quantification];
    B-->E[MultiQC - single HTML report];
    C--Gene-level read counts-->F[DESeq2 differential expression analysis];
    F--normalized counts-->G[ESTIMATE via R package tidyestimate];
    G--adjusted model-->F;
    F--DEG lists-->H[Cluster Profiler];
    F--DEG lists + normalized counts-->I[decoupleR];
    D--TPM counts-->J[Biomarker filtering - custom R scripts];
    F--DEG lists + normalized counts-->J;
    J-->K[ML];
    J--protein-coding genes-->L[STRINGDB];
```
