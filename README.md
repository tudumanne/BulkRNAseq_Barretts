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
    A[Raw FASTQ files]-->B[STAR alignment];
    A-->C[Salmon quantification];
    B-->D[Gene-level read counts];
    D-->E[DESeq2];
    E-->F[Normalized counts];
    E-->G[DEG lists];
    G-->H[Cluster Profiler];
    G-->I[decoupleR]; 
    C-->J[TPM counts];
    G-->K[Biomarker filtering];
    J-->K;
    K-->L[Shortlisted biomarkers]
    L-->M[ML ranking];
    L-->N[protein-coding genes]
    N-->O[STRINGDB];
```
