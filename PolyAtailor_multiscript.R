library(PolyAtailor)
library(tidyverse)
library(Biostrings) 
library(parallel) 
library(dplyr)
library(Rsamtools) 

conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("rbind", "base")
conflict_prefer("cbind", "base")
conflict_prefer("strsplit", "base")
conflict_prefer("count", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("reduce", "IRanges")
conflict_prefer("geom_bar", "ggplot2")
conflict_prefer("first", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("compose", "purrr")
conflict_prefer("last", "dplyr")
conflict_prefer("simplify", "purrr")

#initialize variables using input; input must be 4 digit sample number
#this should work as long as each sample starts its own R session, if bugs, can try tailMap function???
sample_num <- as.character(as.numeric(commandArgs(trailingOnly = TRUE)[1]))

R1name = paste(sample_num, "_1_rc.fastq", sep = "")
R2name = paste(sample_num, "_2.fastq", sep = "")
bamname = paste(sample_num, "_alignmentAligned.sortedByCoord.out.bam", sep = "")
sourcefastq = "reads/"
sourcebam = "alignments/"
resultfolder = "results/"
setwd(resultfolder)

#tailScan parallelized helper function
parallel_helper = function(filename) {
  tailScan(filename,mcans=5,findUmi = F,resultpath = "./",samplename = filename,tailAnchorLen=8,minTailLen=8,realTailLen=20,maxNtail=1,mapping=F)
}

#performs tailScan in 1 direction
perform_tailscan <- function(samplename) {
  
  #create a shortened name for binfiles and output csv
  shortname = substr(samplename, 1, nchar(samplename) - 6)
  binname = paste(shortname, "bin", sep ="")
  
  #unzip and split
  system(paste("zcat ", sourcefastq, samplename, ".gz > ", resultfolder, samplename, sep = ""))
  system(paste("split -l 600000 ", samplename, " ", binname, sep = ""))  
  
  #perform tailScan and prep csv
  filelist <- list.files(pattern = binname)
  output<-mclapply(filelist, parallel_helper, mc.cores = 4)
  
  fixed_df <- Reduce(full_join, output)
  fixed_df <- subset(fixed_df, select = -c(read_type, rt, sample))

  #cleanup and save
  binfiles <- list.files(pattern = binname)
  file.remove(binfiles)
  write.csv(fixed_df, paste(resultfolder, shortname, "_tscan.csv", sep = ""), row.names = FALSE)
  
  system(paste("rm ", samplename, sep = ""))
}

perform_tailscan(R1name)
perform_tailscan(R2name)

tailMap_result <- tailMap(paste(sourcebam, bamname, sep = ""),mcans=5,minTailLen=8,findUmi = F,longRead=F)
write.csv(tailMap_result, paste(resultfolder, sample_num, "_tmap.csv", sep = ""))

