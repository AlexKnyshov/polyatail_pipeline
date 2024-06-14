#!/bin/bash -l

#$ -P casa              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job               
#$ -N pat_stpa      # Give job a name
#$ -pe omp 12
#$ -tc 10

# USAGE: Perform corrections by coroborating PolyAtailor and samToPolyA outputs on all samples indicated in "sample" file. Steps are split into three R scripts, one python script:

module load R/4.1.2
module load samtools
module load bedtools
output_dir=PolyAtailor_corrections

# all the samples IDs to run this script with
sample=$(sed -n ${SGE_TASK_ID}p samples.txt)

# pat_stpa_preprocess.R : extracts a list of reads from polyAtailor. 
Rscript pat_stpa_preprocess.R $sample

# BAM alignment files are filtered to only retain the reads found in PolyAtailor (pat_reads_${sample}.txt) 
samtools view -N ${output_dir}/pat_reads_${sample}.txt -o ${output_dir}/stp1_${sample}_1.bam data/${sample}_1_alignmentAligned.sortedByCoord.out.bam

# Multimapped reads are removed
samtools view -h -q 255 -o ${output_dir}/stp2_${sample}_1.bam ${output_dir}/stp1_${sample}_1.bam

# Filtered BAM file is indexed
samtools index ${output_dir}/stp2_${sample}_1.bam

# bam_length_sclipped.py : given a BAM, finds the length of each read and number of soft clipped bases
python bam_length_sclipped.py ${output_dir}/stp2_${sample}_1.bam

samtools view -N ${output_dir}/pat_reads_${sample}.txt -o ${output_dir}/stp1_${sample}_2.bam data/${sample}_2_alignmentAligned.sortedByCoord.out.bam
samtools view -h -q 255 -o ${output_dir}/stp2_${sample}_2.bam ${output_dir}/stp1_${sample}_2.bam
samtools index ${output_dir}/stp2_${sample}_2.bam
python bam_length_sclipped.py ${output_dir}/stp2_${sample}_2.bam

# pat_stpa_combine.R : performs filtering/coroborates results from samToPolyA and PolyAtailor. Corrects for mismapping of reads with 
# less than 15 non-tail bases. Creates tscan_1_dist_${sample}.bed, tscan_2_dist_${sample}.bed, tscan_both_dist_${sample}.bed, which are used to calculate distances to nearest terminal exon
Rscript pat_stpa_combine.R $sample 

# sort BED files output by pat_stpa_combine.R before using bedtools closest operation
sortBed -i ${output_dir}/tscan_1_dist_${sample}.bed > ${output_dir}/tscan_1_dist_${sample}.sorted.bed
sortBed -i ${output_dir}/tscan_2_dist_${sample}.bed > ${output_dir}/tscan_2_dist_${sample}.sorted.bed
sortBed -i ${output_dir}/tscan_both_dist_${sample}.bed > ${output_dir}/tscan_both_dist_${sample}.sorted.bed

# for each line in each BED file (which corresponds to a read with a poly(A) tail detected), use bedtools closest to find its nearest terminal exon in the reference genome 
# exons_collapse_updated.sorted.bed is the reference genome file. Exons were extracted, ranges were collapsed, the file was 
# converted to a BED file, and then sorted using a similar sortBed command as above

bedtools closest -D a -d -id -a ${output_dir}/tscan_1_dist_${sample}.sorted.bed -b ${output_dir}/exons_collapse_updated.sorted.bed -s > ${output_dir}/tscan_1_calc_${sample}.bed
bedtools closest -D a -d -id -a ${output_dir}/tscan_2_dist_${sample}.sorted.bed -b ${output_dir}/exons_collapse_updated.sorted.bed -s > ${output_dir}/tscan_2_calc_${sample}.bed
bedtools closest -D a -d -id -a ${output_dir}/tscan_both_dist_${sample}.sorted.bed -b ${output_dir}/exons_collapse_updated.sorted.bed -s > ${output_dir}/tscan_both_calc_${sample}.bed

# find_distances.R using the bed files output by bedtools closest, 
Rscript find_distances.R $sample 

# cleanup
rm ${output_dir}/stp1_${sample}_1.bam
rm ${output_dir}/stp1_${sample}_2.bam
rm ${output_dir}/stp2_${sample}_1.bam
rm ${output_dir}/stp2_${sample}_2.bam
rm ${output_dir}/stp2_${sample}_1.bed
rm ${output_dir}/stp2_${sample}_2.bed
rm ${output_dir}/stp2_${sample}_1.bam.bai
rm ${output_dir}/stp2_${sample}_2.bam.bai
rm ${output_dir}/clipped_bases_${sample}_1.csv
rm ${output_dir}/clipped_bases_${sample}_2.csv
rm ${output_dir}/tscan_1_dist_${sample}.bed
rm ${output_dir}/tscan_2_dist_${sample}.bed
rm ${output_dir}/tscan_both_dist_${sample}.bed
rm ${output_dir}/tscan_1_dist_${sample}.sorted.bed
rm ${output_dir}/tscan_2_dist_${sample}.sorted.bed
rm ${output_dir}/tscan_both_dist_${sample}.sorted.bed

rm -f ${output_dir}/*calc*


