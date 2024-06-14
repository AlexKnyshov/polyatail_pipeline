#!/bin/bash -l

#$ -P casa              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job
#$ -N samToPolyA   # Give job a name
#$ -pe omp 4
#$ -tc 30

# USAGE: Run samToPolyA on all samples in sample_list file. Outputs are directed to outputs/

module load samtools
sample_list=samples.txt
samplename=$(sed -n ${SGE_TASK_ID}p ${sample_list})

samtools view ${samplename}_1_alignmentAligned.sortedByCoord.out.bam | samToPolyA.pl --minClipped=8 --minAcontent=0.8 --discardInternallyPrimed --genomeFasta=ensembl_ref.fa - > outputs/${samplename}_1_samtopolya.bed

samtools view ${samplename}_2_alignmentAligned.sortedByCoord.out.bam | samToPolyA.pl --minClipped=8 --minAcontent=0.8 --discardInternallyPrimed --genomeFasta=ensembl_ref.fa - > outputs/${samplename}_2_samtopolya.bed
