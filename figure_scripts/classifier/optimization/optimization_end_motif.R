library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
path <- "data/end_motifs"
healthy_path <- "hbc/end_motifs"
outdir <- ""
title <- "end_motif"
title2 <- "End Motif"

### Import data 
data <- read.delim(list.files(path, "CHARM_LFS_genome2_end_motifs.txt", full.names = TRUE))
normal <- read.delim(list.files(healthy_path, "CHARM_HBC_genome2_end_motifs.txt", full.names = TRUE))
data_samples <- read.delim("sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
data[is.na(data)] <- 0
row.names(data) <- data$motif
data <- data[, colnames(data) %in% data_samples$sWGS]
data <- as.data.frame(t(data))
data <- data[data_samples$sWGS, ]

normal[is.na(normal)] <- 0
row.names(normal) <- normal$motif
normal <- normal[, -1]
normal <- as.data.frame(t(normal))

### Run classifiers
source("figure_scripts/classifier/optimization/runner.R")

