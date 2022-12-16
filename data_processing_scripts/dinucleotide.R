library(tidyverse)
library(dplyr)

### Set working variables
analysis <- "TP53"
path <- file.path("/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/dinucleotide/output", analysis)
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

data_raw <- data.frame(matrix(ncol = 268, nrow = 0))
colnames(data_raw) <- c(colnames(datalist[[1]]), "sample")

for (name in names) {
  x <- datalist[[name]]
  x$sample <- name
  data_raw <- rbind(data_raw, x)
}

### Read in data and combine (Percentages)
datalist <- lapply(filenames, function(x) {read.delim(file = x)})
names(datalist) <- names
data <- data.frame(matrix(ncol = 268, nrow = 0))
colnames(data) <- c(colnames(datalist[[1]]), "sample")

for (name in names) {
  x <- datalist[[name]]
  x$sample <- name
  data <- rbind(data, x)
}

### Write files
write.table(data_raw, file.path(outdir, paste0(project, "_", analysis, "_dinucleotide_raw.txt")), row.names = FALSE, sep = "\t")
write.table(data, file.path(outdir, paste0(project, "_", analysis, "_dinucleotide.txt")), row.names = FALSE, sep = "\t")
