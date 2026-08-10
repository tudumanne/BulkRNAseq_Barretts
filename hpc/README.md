# HPC Pipeline for Bulk RNA-seq Processing

This directory contains the High Performance Computing (HPC) pipeline used to preprocess bulk RNA-seq data prior to downstream analyses in R. The pipeline was run on Linux using Slurm.

## Workflow

``text 
Raw FASTQ files       
      ├──────────────────────├──────────────────────────├
      ▼                      ▼                          ▼
    FastQC             STAR alignment          Salmon quantification
      │                      │                          │ 
      ▼                      ▼                          ▼
    MultiQC        Gene-level read counts          TPM counts
      │
      ▼  
   QC reports 
``

## Software

The following software modules were used for the analysis.

| Software | Version |
|----------|---------|
| FastQC | 0.12.1 |
| MultiQC | 1.13 |
| STAR | 2.7.11b |
| Salmon | 1.10.1 |
