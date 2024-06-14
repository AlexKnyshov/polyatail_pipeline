library(tidyverse)
library(GenomicRanges)
setwd("PolyA")

sample_num <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

corrected_tscan <- read.table(paste("PolyAtailor_corrections/filtered_tscan_", sample_num, ".csv", sep = ""))

# take the filtered tscan results, split into data that has R1 results, R2 results, or R1 and R2 (both) results
coords_both <- corrected_tscan[!is.na(corrected_tscan$start_1) & !is.na(corrected_tscan$start_2),]
coords_1 <- corrected_tscan[!is.na(corrected_tscan$start_1) & is.na(corrected_tscan$start_2),]
coords_2 <- corrected_tscan[is.na(corrected_tscan$start_1) & !is.na(corrected_tscan$start_2),]

# load in data from running BEDtools closest function to calculate distance between nearest terminal exon
# and polyadenylated reads
distance_1 <- read.table(paste("PolyAtailor_corrections/tscan_1_calc_", sample_num, ".bed", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
distance_1 <- data.frame(
  readname = distance_1$V4,
  tscan_chr = distance_1$V1,
  tscan_start = distance_1$V2,
  tscan_end = distance_1$V3,
  x = distance_1$V4,
  tscan_strand = distance_1$V6,
  exon_chr = distance_1$V7,
  exon_start = distance_1$V8,
  exon_end = distance_1$V9,
  gene_id = distance_1$V10,
  z = distance_1$V11,
  exon_strand = distance_1$V12,
  distance = -distance_1$V13
)

distance_2 <- read.table(paste("PolyAtailor_corrections/tscan_2_calc_", sample_num, ".bed", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
distance_2 <- data.frame(
  readname = distance_2$V4,
  tscan_chr = distance_2$V1,
  tscan_start = distance_2$V2,
  tscan_end = distance_2$V3,
  x = distance_2$V4,
  tscan_strand = distance_2$V6,
  exon_chr = distance_2$V7,
  exon_start = distance_2$V8,
  exon_end = distance_2$V9,
  gene_id = distance_2$V10,
  z = distance_2$V11,
  exon_strand = distance_2$V12,
  distance = -distance_2$V13
)

distance_both <- read.table(paste("PolyAtailor_corrections/tscan_both_calc_", sample_num, ".bed", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
distance_both <- data.frame(
  readname = distance_both$V4,
  tscan_chr = distance_both$V1,
  tscan_start = distance_both$V2,
  tscan_end = distance_both$V3,
  x = distance_both$V4,
  tscan_strand = distance_both$V6,
  exon_chr = distance_both$V7,
  exon_start = distance_both$V8,
  exon_end = distance_both$V9,
  gene_id = distance_both$V10,
  z = distance_both$V11,
  exon_strand = distance_both$V12,
  distance = -distance_both$V13
)

distance_1 <- subset(distance_1, select=-c(x, z))
distance_2 <- subset(distance_2, select=-c(x, z))
distance_both <- subset(distance_both, select=-c(x, z))

distances <- rbind(distance_1, distance_2, distance_both)

#remove multi-genes
distances <- distances[!grepl(",", distances$gene_id), ]
duplicates <- duplicated(distances[c("gene_id", "readname")])

singles <- distances[!duplicates,]
singles <- singles[singles$gene_id != ".", ]

write.csv(singles, paste("PolyAtailor_corrections/distances_", sample_num, ".csv", sep = ""), row.names = F)
