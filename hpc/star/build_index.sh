#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --time=24:00:00
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=4
#SBATCH --account=X
#SBATCH --output=../logs/star_%j.out

module load STAR

source ../config.sh

STAR --runMode genomeGenerate \
--genomeDir ${STAR_INDEX} \
--genomeFastaFiles ${GENOME} \
--sjdbOverhang 100 \
--sjdbGTFfile ${GTF} \
--runThreadN ${THREADS}