library(tidyverse)
library(plyr)
library(matrixStats)
library(reshape2)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/hereditary/projects/CHARM/LFS/repeat_masker/output"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/repeat_masker"
project <- "TGL49_LFS"

### Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### Find samples
samples <- list.dirs(path, full.names = FALSE)
samples <- samples[!(samples == "")]

### Find analyses
names <- list.files(file.path(path, samples[[1]]), pattern = "*.txt", full.names = FALSE)
names <- gsub(samples[[1]], "", names)
names <- gsub("_rm_", "", names)
names <- gsub(".txt", "", names)

### Make dataframes
cohort_count <- data.frame(length = rep(10:600, 47),
                           type = rep(names, each = 591))
cohort_freq <- data.frame(length = rep(10:600, 47),
                          type = rep(names, each = 591))
cohort_prop <- data.frame(type = names)

### Loop over each sample
for (i in c(1:length(samples))) {
  sample <- samples[[i]]
  analyses <- list.files(file.path(path, sample), pattern = "*.txt", full.names = TRUE)
  
  ### Read in data
  datalist <- lapply(analyses, function(x){read.delim(file = x, skip = 13, colClasses = c(rep("integer", 2), rep("NULL", 2)), 
                                                      header = TRUE, col.names = c("length", "freq", "X", "Y"))})
  data <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)
  colnames(data) <- c("length", names)
  
  ### Calculate frequency and proportion
  data[is.na(data)] <- 0
  data <- data[data$length %in% c(10:600),]
  
  sums <- colSums(data[, 2:ncol(data)])
  freq <- sweep(data[, 2:ncol(data)], 2, sums, "/")
  freq$length <- data$length
  
  short <- data[data$length %in% c(10:150), ]
  prop <- data.frame(type = names,
                     prop = colSums(short[, 2:ncol(short)])/sums)
  
  ### Merge output to cohort
  data_melt <- reshape2::melt(data, id = "length")
  colnames(data_melt) <- c("length", "type", sample)
  cohort_count <- merge(cohort_count, data_melt, by = c("length", "type"), all = TRUE)
  
  freq_melt <- reshape2::melt(freq, id = "length")
  colnames(freq_melt) <- c("length", "type", sample)
  cohort_freq <- merge(cohort_freq, freq_melt, by = c("length", "type"), all = TRUE)
  
  colnames(prop) <- c("type", sample)
  cohort_prop <- merge(cohort_prop, prop, by = "type", all = TRUE)
}

### Save data
write.table(cohort_count, file.path(outdir, paste0(project, "_rm_counts.txt")), row.names = FALSE, sep = "\t")
write.table(cohort_freq, file.path(outdir, paste0(project, "_rm_freq.txt")), row.names = FALSE, sep = "\t")
write.table(cohort_prop, file.path(outdir, paste0(project, "_rm_prop.txt")), row.names = FALSE, sep = "\t")




