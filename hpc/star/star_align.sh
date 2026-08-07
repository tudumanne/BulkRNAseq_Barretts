#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --time=48:00:00
#SBATCH --mem=40GB
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --account=X
#SBATCH --output=../logs/star_%j.out

module load STAR

source ../config.sh

for i in ${FASTQ}/*_R1.fastq.gz;
  do name=$(basename ${i} _R1.fastq.gz); 
    STAR --readFilesIn ${FASTQ}/${name}_R1.fastq.gz ${FASTQ}/${name}_R2.fastq.gz \
    --genomeDir ${STAR_INDEX} --readFilesCommand zcat --outFileNamePrefix ${STAR}/${name} --runThreadN ${THREADS} --twopassMode Basic \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --outFilterScoreMinOverLread 0.33 \
    --outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 --outSAMstrandField intronMotif --outFilterIntronMotifs None \
    --alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted --outSAMunmapped Within \
    --genomeLoad NoSharedMemory --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
    --chimOutJunctionFormat 1 --chimMainSegmentMultNmax 1 --outSAMattributes NH HI AS nM NM ch;
  done