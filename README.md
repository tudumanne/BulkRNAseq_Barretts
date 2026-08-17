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
    D-->E[DESeq2 (R)];
    E-->F[normalized counts];
    F-->G[tidyestimate (R)];
    G-->H[stromal/immune scores];
    H-->E;
    E-->I[DEG lists]
    I-->J[Cluster Profiler (R)];
    F-->K[decoupleR (R)]; 
    I-->K;
    C-->L[TPM counts];
    E-->M[Biomarker filtering (R)];
    F-->M;
    L-->M;
    M-->K[ML ranking (R)];
    M-->N[protein-coding genes]
    N-->O[STRINGDB];
```
