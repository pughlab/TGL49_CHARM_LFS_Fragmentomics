library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/breakpoint"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/breakpoint"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization"
title <- "breakpoint"
title2 <- "Breakpoint"

### Import data 
data <- read.delim(list.files(path, "CHARM_LFS_genome_breakpoint_ratio.txt", full.names = TRUE))
normal <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_breakpoint_ratio.txt", full.names = TRUE))
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data$sample, ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
data[is.na(data)] <- 0
data <- reshape2::melt(data, id = c("nucleotide", "sample"))
data$context <- paste0(data$nucleotide, "_", data$variable)
data <- reshape2::dcast(data[, c("sample", "value", "context")], context ~ sample)

row.names(data) <- data$context
data <- data[, colnames(data) %in% data_samples$sWGS]
data <- as.data.frame(t(data))
data <- data[data_samples$sWGS, ]

normal[is.na(normal)] <- 0
normal <- reshape2::melt(normal, id = c("nucleotide", "sample"))
normal$context <- paste0(normal$nucleotide, "_", normal$variable)
normal <- reshape2::dcast(normal[, c("sample", "value", "context")], context ~ sample)

row.names(normal) <- normal$context
normal <- normal[, -1]
normal <- as.data.frame(t(normal))

### Run classifiers
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/runner.R")

