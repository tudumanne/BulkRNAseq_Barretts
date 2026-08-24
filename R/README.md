## Downstream analyses using R

This directory contains the R scripts used for the downstream analyses. 


### R packages

The following R Bioconductor packages were used for the analysis.

| Software | Version |
|----------|---------|
| FastQC | 0.12.1 |
| MultiQC | 1.13 |
| STAR | 2.7.11b |
| Salmon | 1.10.1 |

In addition, 


### Input files and scripts

Required inputs include:

- Paired-end FASTQ files (`*.fastq.gz`)
- Reference genome (FASTA, GRCh38.primary_assembly.genome.fa)
- Gene annotation (GTF, gencode.v36.primary_assembly.annotation.gtf)

The scripts required for the analysis are located in respective directories. 
Project-specific paths are defined in

```text
config.sh
```


### Workflow summary

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
