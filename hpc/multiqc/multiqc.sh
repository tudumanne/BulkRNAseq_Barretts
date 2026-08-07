#!/bin/bash
#SBATCH --job-name=MultiQC
#SBATCH --time=2:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --account=X
#SBATCH --output=../logs/multiqc_%j.out

Module load MultiQC

source ../config.sh

mkdir -p ${MULTIQC}

multiqc ${FASTQC} \
-o ${MULTIQC} 