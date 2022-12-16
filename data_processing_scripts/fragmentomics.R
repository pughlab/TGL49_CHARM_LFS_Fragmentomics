library(tidyverse)
library(plyr)
library(vroom)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/fragment_ratio/output"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/fragmentomics"
project <- "CHARM_LFS"

filenames_100kb <- list.files(path = path, pattern = "100kb_bins.txt", recursive = TRUE, full.names = TRUE)
filenames_5mb <- list.files(path = path, pattern = "5Mb_bins.txt", recursive = TRUE, full.names = TRUE)
names <- list.files(path = path, pattern = "5Mb_bins.txt", recursive = TRUE, full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in data
samples <- read.delim(samples)

rows_100kb <- filenames_100kb[1]
rows_100kb <- vroom(rows_100kb)[ , 1:3]

rows_5mb <- filenames_5mb[1]
rows_5mb <- read.delim(rows_5mb, header = TRUE)[ , 2:5]

datalist <- lapply(filenames_100kb, function(x){vroom(file = x, col_select = c(ratio.corrected))})
ratio_100kb <- do.call(cbind, datalist)
ratio_100kb <- as.data.frame(ratio_100kb)

datalist <- lapply(filenames_100kb, function(x){vroom(file = x, col_select = c(nfrags.corrected))})
frags_100kb <- do.call(cbind, datalist)
frags_100kb <- as.data.frame(frags_100kb)

datalist <- lapply(filenames_5mb, function(x){vroom(file = x, col_select = c(ratio_centered))})
ratio_5mb <- do.call(cbind, datalist)
ratio_5mb <- as.data.frame(ratio_5mb)

datalist <- lapply(filenames_5mb, function(x){vroom(file = x, col_select = c(coverage_centered))})
cov_5mb <- do.call(cbind, datalist)
cov_5mb <- as.data.frame(cov_5mb)

### Format fragmentation ratios and write data
colnames(ratio_100kb) <- names
ratio_100kb <- ratio_100kb[ , colnames(ratio_100kb) %in% samples$sWGS]
ratio_100kb <- cbind(rows_100kb, ratio_100kb)
write.table(ratio_100kb, file.path(outdir, paste0(project, "_ratio_100kb.txt")), row.names = FALSE, sep = "\t")

colnames(frags_100kb) <- names
frags_100kb <- frags_100kb[ , colnames(frags_100kb) %in% samples$sWGS]
frags_100kb <- cbind(rows_100kb, frags_100kb)
write.table(frags_100kb, file.path(outdir, paste0(project, "_frags_100kb.txt")), row.names = FALSE, sep = "\t")

colnames(ratio_5mb) <- names
ratio_5mb <- ratio_5mb[ , colnames(ratio_5mb) %in% samples$sWGS]
ratio_5mb <- cbind(rows_5mb, ratio_5mb)
write.table(ratio_5mb, file.path(outdir, paste0(project, "_ratio_5Mb.txt")), row.names = FALSE, sep = "\t")

colnames(cov_5mb) <- names
cov_5mb <- cov_5mb[ , colnames(cov_5mb) %in% samples$sWGS]
cov_5mb <- cbind(rows_5mb, cov_5mb)
write.table(cov_5mb, file.path(outdir, paste0(project, "_coverage_5Mb.txt")), row.names = FALSE, sep = "\t")
