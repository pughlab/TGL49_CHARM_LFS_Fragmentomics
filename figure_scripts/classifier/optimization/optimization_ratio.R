library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization"
title <- "ratio"
title2 <- "Fragment Ratio"

### Import data (Starting with the 5Mb ratios)
data_samples <- read.delim(file.path(path, "samples/sample_list.txt"))
data <- read.delim(file.path(path, "fragmentomics", "CHARM_LFS_ratio_5Mb.txt"))
normal <- read.delim(file.path(healthy, "fragmentomics", "TGL49_HBC_ratio_5Mb.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
row.names(data) <- paste0(data$seqnames, "_", data$start)
data <- data[, data_samples$sWGS]
data <- as.data.frame(t(data))

row.names(normal) <- paste0(normal$seqnames, "_", normal$start)
normal <- normal[, !(colnames(normal) %in% c("seqnames", "arm", "start", "end"))]
normal <- as.data.frame(t(normal))

### Run classifiers
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/runner.R")
