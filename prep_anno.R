library(tidyverse)
library(GenomicRanges)
library("rtracklayer")
gtf_data = import('reference/Homo_sapiens.GRCh38.108.gtf')

exons <- gtf_data[gtf_data$type == "exon"]

# split exons by transcript ID
exons_by_transcript <- split(exons, mcols(exons)$transcript_id)

# create GRanges object to store last exons
last_exons <- GRanges()

# iterate through each transcript
for (transcript_id in names(exons_by_transcript)) {
  # create GRanges object for each transcript, sort by coordinates 
  sample <- exons_by_transcript[[transcript_id]]
  
  if (runValue(strand(sample[1]) == "+")) {
    transcript_exons <- sample[order(start(sample))]
    }
  else {
    transcript_exons <-sample[order(start(sample), decreasing = T)]
    }
  
  # find the last exon for the current transcript
  cur_exon <- transcript_exons[length(transcript_exons)]
  
  # append the last exon to the 'last_exons' object
  last_exons <- append(last_exons, cur_exon)
}

#recover regions flanking 10k bp downstream the set of ranges in last_exons
#final <- resize(last_exons, width=width(last_exons) + 10000)

# Write the last exons to a new GTF file
export(last_exons, "exons_noflank_updated.gtf")

