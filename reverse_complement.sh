#!/bin/bash -l

#$ -P casa              # Specify the SCC project name you want to use
#$ -l h_rt=2:00:00     # Specify the hard time limit for the job               
#$ -N rev_comp          # Give job a name
#$ -pe omp 12

# USAGE: Generate reverse complement versions of R1 fastq files for all samples indicated in the "sample" file. Outputs directed to data/

# setup array job
sample=$(sed -n ${SGE_TASK_ID}p samples.txt)

zcat location/${sample}_1.fastq.gz | seqkit seq -r -p -t RNA | gzip > data/${sample}_1_rc.fastq.gz
