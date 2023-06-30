library(tidyverse)
library(dplyr)

### Set working variables
path <- file.path("/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/dinucleotide/output/short")
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/dinucleotide"
project <- "CHARM_LFS"

filenames <- list.files(path = path, pattern = "contexts.txt", recursive = TRUE, full.names = TRUE)
filenames_raw <- list.files(path = path, pattern = "raw.txt", recursive = TRUE, full.names = TRUE)
names <- list.files(path = path, pattern = "contexts.txt", recursive = TRUE, full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in data and combine (Raw counts)
datalist <- lapply(filenames_raw, function(x) {read.delim(file = x)})
names(datalist) <- names

data_raw <- Reduce(function(x,y){merge(x,y, by = "context", all = TRUE)}, datalist)
colnames(data_raw) <- c("context", names)

### Read in data and combine (Percentages)
datalist <- lapply(filenames, function(x) {read.delim(file = x)})
names(datalist) <- names

data <- Reduce(function(x,y){merge(x,y, by = "context", all = TRUE)}, datalist)
colnames(data) <- c("context", names)

### Write files
write.table(data_raw, file.path(outdir, paste0(project, "_short_dinucleotide_raw.txt")), row.names = FALSE, sep = "\t")
write.table(data, file.path(outdir, paste0(project, "_short_dinucleotide.txt")), row.names = FALSE, sep = "\t")
