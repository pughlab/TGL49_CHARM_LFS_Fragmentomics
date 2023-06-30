library(tidyverse)
library(dplyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/hereditary/projects/CHARM/LFS/breakpoint/output"
outdir <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/breakpoint"
project <- "CHARM_LFS"
analyses <- list.dirs(path, full.names = FALSE)
analyses <- analyses[!(analyses == "")]

for (analysis in analyses) {
  site_path <- file.path(path, analysis)
  filenames <- list.files(path = site_path, pattern = "count.txt", recursive = TRUE, full.names = TRUE)
  names <- list.files(path = site_path, pattern = "count.txt", recursive = TRUE, full.names = FALSE)
  names <- sub("(_WG).*", '\\1', names)
  
  ### Make outdir
  dir.create(outdir, showWarnings = FALSE)
  
  ### Read in data
  datalist <- lapply(filenames, function(x) {read.delim(file = x)})
  names(datalist) <- names
  
  ### Process and combine
  data_count <- data.frame()
  data_freq <- data.frame()
  data_ratio <- data.frame()
  
  for (i in 1:length(datalist)) {
    name <- names[[i]]
    
    x <- datalist[[i]]
    row.names(x) <- x$nucleotide
    
    sum <- sum(x[, 2])
    y <- x[, -1]/sum
    
    z <- y/c(0.2951856, 0.2039184, 0.2047607, 0.2961353) ### Frequencies of each nucleotide genome-wide
    
    x$sample <- name
    data_count <- rbind(data_count, x)
    
    y$nucleotide <- row.names(y)
    y$sample <- name
    data_freq <- rbind(data_freq, y)
    
    z$nucleotide <- row.names(z)
    z$sample <- name
    data_ratio <- rbind(data_ratio, z)
  }
  
  ### Write files
  write.table(data_count, file.path(outdir, paste0(project, "_", analysis, "_breakpoint_count.txt")), row.names = FALSE, sep = "\t")
  write.table(data_freq, file.path(outdir, paste0(project, "_", analysis, "_breakpoint_freq.txt")), row.names = FALSE, sep = "\t")
  write.table(data_ratio, file.path(outdir, paste0(project, "_", analysis, "_breakpoint_ratio.txt")), row.names = FALSE, sep = "\t")
}
