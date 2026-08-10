## HPC Pipeline for Bulk RNA-seq Processing

This directory contains the High Performance Computing (HPC) pipeline used to preprocess bulk RNA-seq data prior to downstream analyses in R. The pipeline was run on Linux using Slurm.

### Software

The following software modules were used for the analysis.

| Software | Version |
|----------|---------|
| FastQC | 0.12.1 |
| MultiQC | 1.13 |
| STAR | 2.7.11b |
| Salmon | 1.10.1 |

### Workflow

```mermaid
flowchart TD 
    A [Raw FASTQ files] --> B [FastQC]
    A [Raw FASTQ files] --> C [STAR alignment]
    A [Raw FASTQ files] --> D [Salmon quantification]
    B --> E [MultiQC]
    E --> F [QC reports]
    C --> G [Gene-level read counts]
    D --> [TPM counts]
```





