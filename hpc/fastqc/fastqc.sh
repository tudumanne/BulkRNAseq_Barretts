#!/bin/bash
#SBATCH --job-name=FastQC
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=4
#SBATCH --account=X
#SBATCH --output=../logs/fastqc_%j.out

module load FastQC

source ../config.sh

mkdir -p ${FASTQC}

fastqc ${FASTQ}/*.fastq.gz \
--outdir ${FASTQC} \
--nogroup \
--format fastq \
-t ${THREADS} 