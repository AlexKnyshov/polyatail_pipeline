library(tidyverse)
setwd("PolyA")

sample_num <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

#load polyAtailor and samToPolyA data 

pat_R1 <- read.csv(paste("data/", sample_num, "_1_rc_tscan.csv", sep=""))
pat_R1 <- pat_R1[-which(pat_R1$polyAT == "-"),]

pat_R2 <- read.csv(paste("data/", sample_num, "_2_tscan.csv", sep=""))
pat_R2 <- pat_R2[-which(pat_R2$polyAT == "-"),]

#extract readnames from polyAtailorData for bam filtering

pat_merged <- merge(pat_R1, pat_R2, by = "read_num", all = T, suffixes = c("_R1","_R2"))
colnames(pat_merged)[colnames(pat_merged) == "read_num"] <- "readname"

pat_reads <- pat_merged$readname
writeLines(pat_reads, file(paste("PolyAtailor_corrections/pat_reads_", sample_num, ".txt", sep = "")))
