#!/bin/bash

########## project directories ##########

PROJECT = /path/to/project

FASTQ = ${PROJECT}/data/raw_fastq
REFERENCE = ${PROJECT}/reference_genome
OUTPUTS = ${PROJECT}/outputs

FASTQC = ${OUTPUTS}/fastqc
MULTIQC = ${OUTPUTS}/multiqc
STAR = ${OUTPUTS}/star
SALMON = ${OUTPUTS}/salmon

########## reference ##########

GENOME = ${REFERENCE}/GRCh38.primary_assembly.genome.fa
GTF = ${REFERENCE}/gencode.v36.primary_assembly.annotation.gtf
STAR_INDEX = ${REFERENCE}/STAR_index
SALMON_INDEX = ${REFERENCE}/SALMON_index
GENE_ID_LIST = ${REFERENCE}/gene_id_list

########## threads ##########

THREADS = 16
