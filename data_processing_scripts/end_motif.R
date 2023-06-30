library(tidyverse)
library(dplyr)

### Set working variables
path <- file.path("/Users/derekwong/Desktop/H4H/hereditary/projects/CHARM/LFS/end_motifs/output")
outdir <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/end_motifs"
project <- "CHARM_LFS"

analyses <- list.dirs(path, full.names = FALSE)
analyses <- analyses[!(analyses == "")]

for (analysis in analyses) {
  site_path <- file.path(path, analysis)
  filenames <- list.files(path = site_path, pattern = "motifs.txt", recursive = TRUE, full.names = TRUE)
  filenames_raw <- list.files(path = site_path, pattern = "raw.txt", recursive = TRUE, full.names = TRUE)
  names <- list.files(path = site_path, pattern = "motifs.txt", recursive = TRUE, full.names = FALSE)
  names <- sub("(_WG).*", '\\1', names)
  
  ### Make outdir
  dir.create(outdir, showWarnings = FALSE)
  
  ### Read in data and combine (Raw counts)
  datalist <- lapply(filenames_raw, function(x) {read.delim(file = x)})
  names(datalist) <- names
  
  data_raw <- Reduce(function(x, y) merge(x, y, by = "motif", all = TRUE), datalist) 
  colnames(data_raw) <- c("motif", names)
  
  ### Read in data and combine (Percentages)
  datalist <- lapply(filenames, function(x) {read.delim(file = x)})
  names(datalist) <- names
  
  data <- Reduce(function(x, y) merge(x, y, by = "motif", all = TRUE), datalist) 
  colnames(data) <- c("motif", names)
  
  ### Write files
  write.table(data_raw, file.path(outdir, paste0(project, "_", analysis, "_end_motifs_raw.txt")), row.names = FALSE, sep = "\t")
  write.table(data, file.path(outdir, paste0(project, "_", analysis, "_end_motifs.txt")), row.names = FALSE, sep = "\t")
}
