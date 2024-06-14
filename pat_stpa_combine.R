library(tidyverse)
library(GenomicRanges)
setwd("PolyA")

sample_num <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

# load in samToPolyA data
stpa_R1 <- read.table(paste("samToPolyA_alldata/", sample_num, "_1_samtopolya.bed", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
stpa_R2 <- read.table(paste("samToPolyA_alldata/", sample_num, "_2_samtopolya.bed", sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# create a dataframe from the BED data
stpa_R1 <- data.frame(
  chromosome = stpa_R1$V1,
  start = stpa_R1$V2,
  end = stpa_R1$V3,
  readname = stpa_R1$V4,
  PAL = stpa_R1$V5,
  strand = stpa_R1$V6
)

stpa_R2 <- data.frame(
  chromosome = stpa_R2$V1,
  start = stpa_R2$V2,
  end = stpa_R2$V3,
  readname = stpa_R2$V4,
  PAL = stpa_R2$V5,
  strand = stpa_R2$V6
)

pat_R1 <- read.csv(paste("PolyAtailor_alldata/", sample_num, "_1_rc_tscan.csv", sep=""))
pat_R1 <- pat_R1[-which(pat_R1$polyAT == "-"),]

pat_R2 <- read.csv(paste("PolyAtailor_alldata/", sample_num, "_2_tscan.csv", sep=""))
pat_R2 <- pat_R2[-which(pat_R2$polyAT == "-"),]

pat_merged <- merge(pat_R1, pat_R2, by = "read_num", all = T, suffixes = c("_R1","_R2"))
colnames(pat_merged)[colnames(pat_merged) == "read_num"] <- "readname"

# function to calculate the length of the longest consecutive 'T' sequence with 1 non-'T' character allowed
calculate_longest_t_sequence <- function(sequence) {
  tolerance = 1
  if (is.na(sequence)) {
    return(NA)
  }
  segments <- unlist(strsplit(sequence, "(?<=T)(?=[^T])|(?<=[^T])(?=T)", perl=TRUE))
  num_segments <- length(segments)
  
  if (num_segments == 1) {
    return(nchar(sequence))
  }
  if (num_segments == 2) {
    return(nchar(segments[1]) + as.integer(min(segments[2], tolerance)))
  }
  
  max_length <- nchar(max(segments))
  total_length <- 0
  
  for (i in seq(1, num_segments-2, by = 2)) {
    if (nchar(segments[i+1]) == tolerance) {
      total_length <- nchar(segments[i]) + nchar(segments[i+1]) + nchar(segments[i+2])
      if (total_length > max_length) {
        max_length <- total_length
      }
    }
  }
  
  return(max_length)
}

# adjust the PAL to account for one mismatch, using the function above
pat_merged$t1_R1 <- sapply(pat_merged$tail_R1, calculate_longest_t_sequence)
pat_merged$t1_R2 <- sapply(pat_merged$tail_R2, calculate_longest_t_sequence)

# load in the bedfile data containing the coordinates of each read, excluding multimapped reads
bed_1 <- read.table(paste("PolyAtailor_corrections/stp2_", sample_num, "_1.bed", sep=""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

bed_1 <- data.frame(
  chromosome = bed_1$V1,
  start = bed_1$V2,
  end = bed_1$V3,
  readname = bed_1$V4,
  strand = bed_1$V5,
  length = bed_1$V6
)

pat_merged_1 <- merge(pat_merged, bed_1, by = "readname", all.x = TRUE)

colnames(pat_merged_1) <- c("readname", "polyAT_R1", "PAL_R1", "tail_R1", "tailTypeR1", "nA_R1", "polyAT_R2", "PAL_R2", "tail_R2", "tailType_R2", "nA_R2", "t1_R1", "t1_R2", "chr_1", "start_1", "end_1", "strand_1", "length_1")

bed_2 <- read.table(paste("PolyAtailor_corrections/stp2_", sample_num, "_2.bed", sep=""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)

bed_2 <- data.frame(
  chromosome = bed_2$V1,
  start = bed_2$V2,
  end = bed_2$V3,
  readname = bed_2$V4,
  score = bed_2$V5,
  strand = bed_2$V6
)

pat_full <- merge(pat_merged_1, bed_2, by = "readname", all.x = TRUE)

colnames(pat_full) <- c("readname", "polyAT_R1", "PAL_R1", "tail_R1", "tailTypeR1", "nA_R1", "polyAT_R2", "PAL_R2", "tail_R2", "tailType_R2", "nA_R2", "t1_R1", "t1_R2", "chr_1", "start_1", "end_1", "strand_1", "length_1", "chr_2", "start_2", "end_2", "strand_2", "length_2" )

#separate tailscan data based on R1 only, R2 only, and both R1/R2 data
empty_rows <- which(is.na(pat_full$start_1) & is.na(pat_full$start_2))
pat_full <- pat_full[-empty_rows, ]

both_tailed <- which(!is.na(pat_full$PAL_R1) & !is.na(pat_full$PAL_R2))
one_tail_tscan <- pat_full[-both_tailed, ]
both_tail_tscan <- pat_full[both_tailed,]

R1_tails <- which(is.na(one_tail_tscan$PAL_R2))

# look at reads that have a tail detected in R1 xor R2 and filter
R1_tscan <- one_tail_tscan[R1_tails,]
R2_tscan <- one_tail_tscan[-R1_tails, ]

# create a column that indicates whether the read is found in samtoPolyA data
R1_tscan$sam <- ifelse(R1_tscan$readname %in% unique(stpa_R1$readname), TRUE, FALSE)
R2_tscan$sam <- ifelse(R2_tscan$readname %in% unique(stpa_R2$readname), TRUE, FALSE)

# remove R1 tails that have no R1 coord data, we can't asses whether the polyA transcript is too long in these
# same for R2
R1_tscan <- R1_tscan[-which(is.na(R1_tscan$length_1)),] 
R2_tscan <- R2_tscan[-which(is.na(R2_tscan$length_2)),]

R1_tscan$keep <- FALSE
R1_tscan$pref_coord <- NA

for (i in 1:nrow(R1_tscan)) {
  # if the R1 PAL is too big, use the R2 coords
  if ((R1_tscan$length_1[i] - R1_tscan$t1_R1[i]) < 15) {
    R1_tscan$start_1[i] <- NA
    R1_tscan$end_1[i] <- NA
    R1_tscan$strand_1[i] <- NA
    R1_tscan$chr_1[i] <- NA
    R1_tscan$keep[i] <- TRUE
    R1_tscan$pref_coord[i] <- "R2"
    #R1_tscan$length_1[i] <- NA
  }
  
  else {
    # if the readname is also in samtopolyA, save R1 coords
    if (R1_tscan$sam[i] == TRUE) {
      #R1_tscan$start_2[i] <- NA
      #R1_tscan$end_2[i] <- NA
      #R1_tscan$strand_2[i] <- NA
      #R1_tscan$chr_2[i] <- NA
      R1_tscan$keep[i] <- TRUE
      R1_tscan$pref_coord <- "R1"
      #R1_tscan$length_2[i] <- NA
    }
  
  }
  
} 

R2_tscan$keep <- FALSE
R2_tscan$pref_coord <- NA

for (i in 1:nrow(R2_tscan)) {
  # if the R2 PAL is too big, use the R1 coords, force NA R2 coords
  if ((R2_tscan$length_2[i] - R2_tscan$t1_R2[i]) < 15) {
    R2_tscan$start_2[i] <- NA
    R2_tscan$end_2[i] <- NA
    R2_tscan$strand_2[i] <- NA
    R2_tscan$chr_2[i] <- NA
    R2_tscan$keep[i] <- TRUE
    R2_tscan$pref_coord[i] <- "R1"
    #R2_tscan$length_2[i] <- NA
  }
  
  else {
    # if the readname is also in samtopolyA, save R2 coords, force NA R1 coords
    if (R2_tscan$sam[i] == TRUE) {
      #R2_tscan$start_1[i] <- NA
      #R2_tscan$end_1[i] <- NA
      #R2_tscan$strand_1[i] <- NA
     # R2_tscan$chr_1[i] <- NA
      R2_tscan$keep[i] <- TRUE
      R2_tscan$pref_coord[i] <- "R2"
      #R2_tscan$length_1[i] <- NA
    }
  }
}


# R1_filterd <- R1_tscan[-which(is.na(R1_tscan$readname)),] 
#R1_filterd <- R1_tscan[!(is.na(R1_tscan$start_1) & is.na(R1_tscan$start_2)), ]
R1_filtered <- R1_tscan[(R1_tscan$keep == TRUE), ]
R1_filtered <- R1_filtered[!(is.na(R1_filtered$start_2) & is.na(R1_filtered$start_1)),]

# R2_filterd <- R2_tscan[-which(is.na(R2_tscan$readname)),] 
#R2_filterd <- R2_tscan[!(is.na(R2_tscan$start_1) & is.na(R2_tscan$start_2)), ]
R2_filtered <- R2_tscan[(R2_tscan$keep == TRUE), ]
R2_filtered <- R2_filtered[!(is.na(R2_filtered$start_2) & is.na(R2_filtered$start_1)),]

both_tail_tscan <- both_tail_tscan[!((both_tail_tscan$length_1 - both_tail_tscan$t1_R1) < 15 & (both_tail_tscan$length_2 - both_tail_tscan$t1_R2 < 15)), ]

# take each readname, and if it is in both samToPolyA R1 and R2, and the strands & chrs match in both samtopolyA directions, 
# keep this read and record the coords of both tails
stpa_merge <- merge(stpa_R1, stpa_R2, by = "readname", all = T, suffixes = c("_R1","_R2"))
stpa_merge <- na.omit(stpa_merge)
stpa_merge <- stpa_merge[(stpa_merge$chromosome_R1 == stpa_merge$chromosome_R2), ]
stpa_merge <- stpa_merge[(stpa_merge$strand_R1 == stpa_merge$strand_R2), ]

both_tail_filtered <- both_tail_tscan[both_tail_tscan$readname %in% stpa_merge$readname, ]
both_tail_filtered$pref_coord <- "both"
both_tail_filtered <- both_tail_filtered[!is.na(both_tail_filtered$start_1) & !is.na(both_tail_filtered$start_2),]

R1_filtered <- subset(R1_filtered, select = -sam)
R2_filtered <- subset(R2_filtered, select = -sam)
R1_filtered <- subset(R1_filtered, select = -keep)
R2_filtered <- subset(R2_filtered, select = -keep)

corrected_tscan <- rbind(both_tail_filtered, R1_filtered, R2_filtered)

# if both reads mapped, check if same chr and opposite strand
corrected_tscan <- corrected_tscan[is.na(corrected_tscan$strand_1) | is.na(corrected_tscan$strand_2) | (!is.na(corrected_tscan$strand_1) & !is.na(corrected_tscan$strand_2) & 
                                     corrected_tscan$strand_1 != corrected_tscan$strand_2 & (corrected_tscan$chr_1 == corrected_tscan$chr_2)),]
# remove MT reads
corrected_tscan <- corrected_tscan[(corrected_tscan$chr_1 != "MT" & !is.na(corrected_tscan$chr_1)| corrected_tscan$chr_2 != "MT" & !is.na(corrected_tscan$chr_2)), ]

# add count of soft clipped bases
clipped_1<- read.csv(paste("PolyAtailor_corrections/clipped_bases_", sample_num, "_1.csv", sep = ""))
clipped_2<- read.csv(paste("PolyAtailor_corrections/clipped_bases_", sample_num, "_2.csv", sep = ""))

combined <- merge(clipped_1, clipped_2, by = "readname", suffixes = c("_R1", "_R2"))
soft_clipped <- merge(corrected_tscan, combined, by = "readname", all.x = TRUE)

write.csv(soft_clipped, paste("PolyAtailor_corrections/filtered_tscan_", sample_num, ".csv", sep = ""))

#####

# separate corrected tscan data based on which set of coords are availiable
coords_both <- corrected_tscan[!is.na(corrected_tscan$start_1) & !is.na(corrected_tscan$start_2),]
coords_1 <- corrected_tscan[!is.na(corrected_tscan$start_1) & is.na(corrected_tscan$start_2),]
coords_2 <- corrected_tscan[is.na(corrected_tscan$start_1) & !is.na(corrected_tscan$start_2),]

# convert tscan data to bed file for distance calculation

coords_1_bed <- subset(coords_1, select=c('chr_1', 'start_1', 'end_1', 'readname', 'strand_1', 'strand_1'))
colnames(coords_1_bed) <- c('chrom', 'chromStart', 'chromEnd', 'readname', 'y', 'strand')
write.table(coords_1_bed, paste("PolyAtailor_corrections/tscan_1_dist_", sample_num, ".bed", sep = ""), row.names = F, sep = '\t', quote = F, col.names = F)

coords_2_bed <- subset(coords_2, select=c('chr_2', 'start_2', 'end_2', 'readname', 'strand_2', 'strand_2'))
colnames(coords_2_bed) <- c('chrom', 'chromStart', 'chromEnd', 'readname', 'y', 'strand')
write.table(coords_2_bed, paste("PolyAtailor_corrections/tscan_2_dist_", sample_num, ".bed", sep = ""), row.names = F, sep = '\t', quote = F, col.names = F)

coords_both_bed <- subset(coords_both, select=c('chr_2', 'start_2', 'end_2', 'readname', 'strand_2', 'strand_2'))
colnames(coords_both_bed) <- c('chrom', 'chromStart', 'chromEnd', 'readname', 'y', 'strand')
write.table(coords_both_bed, paste("PolyAtailor_corrections/tscan_both_dist_", sample_num, ".bed", sep = ""), row.names = F, sep = '\t', quote = F, col.names = F)
