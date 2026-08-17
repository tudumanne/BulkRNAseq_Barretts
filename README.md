## RNA-seq Analysis Pipeline for Pathway Analysis and Biomarker Discovery in Barrett's Neoplasia

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
    classDef input fill:#ffffff,stroke:#000000;
    classDef output1 fill:#ffffff,stroke:#000000;
    classDef output2 fill:#edabab,stroke:#520612;
    
    A[Raw FASTQ files]:::input-->B[STAR alignment]:::method;
    A-->C[Salmon quantification]:::method;
    B-->D[Gene-level read counts]:::output1;
    D-->E[DESeq2]:::method;
    E-->F[DEG lists/Normalized counts]:::output1;
    F-->G[Cluster Profiler]:::method;
    F-->H[decoupleR]:::method; 
    C-->I[TPM counts]:::output1;
    F-->J[Biomarker filtering]:::method;
    I-->J;
    J-->K[Shortlisted biomarkers]:::output1;
    K-->L[ML ranking]:::method;
    K-->M[protein-coding genes]:::output1;
    M-->N[STRINGDB]:::method;
    L-->O[five-gene panel]:::output2;
    L-->P[IHC candidates]:::output2;
```
