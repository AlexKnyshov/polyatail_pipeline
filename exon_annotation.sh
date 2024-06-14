#!/bin/bash -l

#$ -P casa              # Specify the SCC project name you want to use
#$ -l h_rt=16:00:00     # Specify the hard time limit for the job               
#$ -N annotation      # Give job a name

# USAGE: create a subset of the original GTF, consisting only of terminal exons

module load R/4.1.2

Rscript prep_anno.R  