#POLY(A) PIPELINE WORKFLOW
Sasha Bacot, modified by Alexander Knyshov
summer 2023

## Pipeline workflow
list of samples used for array job setup on an SGE cluster: `samples.txt`
```
sample1
sample2
sample3
etc
```


1. **Prep files:** Reverse complement R1 FASTQ, split BAMS into R1 and R2, create terminal exon annotation
  * reverse complement script: `reverse_complement.sh`
  * BAM splitting script: `samtools_split.sh`
  * teriminal exon annotation script: `exon_annotation.sh`, which runs `prep_anno.R`
    * after this, 
      * modify your gtf to only retain gene_id in column 9 `sed 's/gene_id [^"]*"\([^"]*\)".*/\1/' FNAME1 > FNAME2`
      * then take that file and pass it to bedtools sort `bedtools sort -i FNAME`
      * then pass it to bedtools merge `bedtools merge -s -i FNAME -c 9,9,7 -o distinct,count_distinct,distinct`
2. **Run PolyAtailor**
  * script: `PolyAtailor.sh`, which runs `PolyAtailor_multiscript.R`
3. **Run samToPolyA**
  * script: `samToPolyA.sh`
4. **Combine PolyAtailor and samToPolyA, perform corrections**
  * script: `pat_stpa_corrections.sh`, this script runs 3 R scripts and 1 python program:
    + `pat_stpa_preprocess.R` takes the polyAtailor outputs, and creates a .txt list of reads. The BAM file is then filtered to only include the reads in this .txt. The BAM is further filtered to remove multimapped reads.
    + `bam_length_sclipped.py` takes the filtered BAM, finds length of each read and number of soft clipped bases.
    + `pat_stpa_combine.R`: performs filtering/coroborates results from samToPolyA and PolyAtailor. Corrects for mismapping of reads withless than 15 non-tail bases. Creates tscan_1_dist_{sample}.bed, tscan_2_dist_{sample}.bed, tscan_both_dist_{sample}.bed, which are used to calculate distances to nearest terminal exon
    + `find_distances.R`: after running Bedtools closest function to find closest terminal exon of each polyadenylated read, this script removes ranges that contain multiple genes after collapsing the exon annotation
5. **Final output generation**
  * count number of uniqely mapped reads to normalize data to CPM: `count_reads.sh`
  * use the R scripts below to create the matrix of cpm values, the matrix of median distance, and to select a set of high confidence polyadenylated genes:
```
library(tidyverse)
library(dplyr)
sample_numbers <- readLines("samples.txt")
exons <- read.table("PolyAtailor_corrections/exons_collapse_updated.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
exons <- data.frame(
  gene_id = exons$V4
)
genes <- exons[!grepl(",", exons$gene_id), ]
genes <- unique(genes)
count_all_genes <- data.frame(gene_id = genes)
gene_distances <- data.frame(gene_id = genes)
# loop through each sample
for (sample_num in sample_numbers) {
  #read input
  file_path <- paste0("PolyAtailor_corrections/distances_", sample_num, ".csv")
  data <- read.csv(file_path)
  #get counts
  gene_counts <- data %>% group_by(gene_id) %>% summarise(read_count = n_distinct(readname))
  count_all_genes <- left_join(count_all_genes, gene_counts, by = "gene_id")
  colnames(count_all_genes)[length(count_all_genes)] <- sample_num
  #get distances
  distances <- data %>% group_by(gene_id) %>% summarise(dist = median(distance))
  gene_distances <- left_join(gene_distances, distances, by = "gene_id")
  colnames(gene_distances)[length(gene_distances)] <- sample_num
}

# convert results to matrix
row_names <- count_all_genes[, 1]
col_names <- colnames(count_all_genes[,-1])
count_all_genes <- as.matrix(count_all_genes[, -1])
rownames(count_all_genes) <- row_names
colnames(count_all_genes) <- col_names
gene_distances <- as.matrix(gene_distances[, -1])
rownames(gene_distances) <- row_names
colnames(gene_distances) <- col_names

# remove genes with no expression across samples from count and dist matrix
count_exp_genes <- count_all_genes[rowSums(!is.na(count_all_genes)) > 0, ]
count_exp_genes[is.na(count_exp_genes)] <- 0
exp_gene_distances <- gene_distances[rowSums(!is.na(gene_distances)) > 0, ]

#import total read counts for scaling to CPM
counts <- read.csv("/restricted/projectnb/amp-ad/sbacot/PolyA/final_outputs/readcounts.csv")
counts$mil <- (counts$count/1000000)/2
cpm_exp_genes <- count_exp_genes
for (i in 1:length(count_exp_genes[1,])) {
  col_name <- colnames(count_exp_genes)[i]
  cpm_exp_genes[,i] <- count_exp_genes[,i] / counts$mil[match(col_name, counts$sample)]
}

polyAsiteHighConfGenes <- rownames(cpm_exp_genes)[apply(cpm_exp_genes, 1, function (x) {sum(x>=1/35)>=128*0.75}) & apply(exp_gene_distances, 1, function (x) {median(x, na.rm=T)==0})]
```
  

