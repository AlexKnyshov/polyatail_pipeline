#!/bin/bash -l

#$ -P casa              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job
#$ -N samtools_split           # Give job a name

# USAGE: Split BAM alignment files into separate R1 and R2 files

module load samtools
samtools view -b -f 0x0040 location/sample_Aligned.sortedByCoord.out.bam > data/sample_R1_PA.bam
samtools index data/sample_R1_PA.bam

samtools view -b -f 0x0080 location/sample_Aligned.sortedByCoord.out.bam > data/sample_R2_PA.bam
samtools index data/sample_R2_PA.bam
