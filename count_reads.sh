#!/bin/bash -l

#$ -P casa              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job               
#$ -N count_reads      # Give job a name
#$ -pe omp 12

# USAGE: counts the number of uniquely mapped reads in every sample indicated by the "sample" file. Outputs results to the "output" csv

module load samtools
output=readcounts.csv

sample=$(sed -n ${SGE_TASK_ID}p samples.txt)

bam_file=data/${sample}_alignmentAligned.sortedByCoord.out.bam
result=$(samtools view -c -q 255 $bam_file)
echo "$sample,$result" >> $output
