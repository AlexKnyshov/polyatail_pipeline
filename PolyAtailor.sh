#!/bin/bash -l

#$ -P casa              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job               
#$ -N PolyAtailor      # Give job a name
#$ -pe omp 12
#$ -tc 10

# USAGE: runs PolyAtailor on all FASTQs (reverse complement R1, normal R2)

module load R/4.1.2

# the maximum number of the sample to run
MAX_SAMPLE=$(sed -n ${SGE_TASK_ID}p samples.txt)

#run R1 and R2 tailScan on the specified sample
Rscript PolyAtailor_multiscript.R $MAX_SAMPLE 
