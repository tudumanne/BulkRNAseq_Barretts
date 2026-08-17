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
    classDef method fill:#d4edda,stroke:#28a745;
    classDef input fill:#ffffff,stroke:#ffffff;
    classDef output1 fill:#ffffff,stroke:#ffffff;
    classDef output2 fill:#edabab,stroke:#edabab;
    
    A[Raw FASTQ files]:::input-->B[STAR alignment]:::method;
    A-->C[Salmon quantification];
    B-->D[Gene-level read counts];
    D-->E[DESeq2];
    E-->F[DEG lists/Normalized counts];
    F-->G[Cluster Profiler];
    F-->H[decoupleR]; 
    C-->I[TPM counts];
    F-->J[Biomarker filtering];
    I-->J;
    J-->K[Shortlisted biomarkers]
    K-->L[ML ranking];
    K-->M[protein-coding genes]
    M-->N[STRINGDB];
    L-->O[five-gene panel];
    L-->P[IHC candidates];
```
