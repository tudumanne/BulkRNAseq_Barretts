## Downstream analyses using R

This directory contains the R scripts used for the downstream analyses. 

### R packages

The following R Bioconductor packages were used for the analyses.

| Package | Version |
|----------|---------|
|  |  |
|  |  |
|  |  |
|  |  |

In addition, 


### Input files and scripts

Required inputs include:

- 
- 
- Gene annotation ()

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
    
    A[Gene-level count data]:::input-->B[DESeq2]:::method;
    B-->C[DEG lists]:::output1;
    B-->D[Normalized counts]:::output1;
    C-->E[Cluster Profiler]:::method;
    C-->F[decoupleR]:::method;
    G[TPM counts]:::input-->H[Biomarker filtering]:::method;
    D-->H;
    C-->H;
    H-->J[Shortlisted biomarkers]:::output1;
    J-->K[ML ranking]:::method;
    J-->L[protein-coding genes]:::output1;
    K-->M[five-gene panel]:::output2;
    K-->N[IHC candidates]:::output2;
```
