#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --time=12:00:00
#SBATCH --mem=40GB
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --account=X
#SBATCH --output=../logs/salmon_%j.out

module load STAR

source ../config.sh

for i in ${FASTQ}/*_R1.fastq.gz;
  do name=$(basename ${i} _R1.fastq.gz); 
    salmon quant -i ${SALMON_INDEX} -l A -1 ${FASTQ}/${name}_R1.fastq.gz -2 ${FASTQ}/${name}_R2.fastq.gz \
    -p ${THREADS} --validateMappings -g ${GENE_ID_LIST} -o ${SALMON}/${name}_quant;
  done