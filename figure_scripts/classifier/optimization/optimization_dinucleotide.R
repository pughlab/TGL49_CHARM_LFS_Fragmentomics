library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(data.table)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/dinucleotide"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/dinucleotide"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization"
title <- "dinucleotide"
title2 <- "Dinucleotide"

### Import data (Starting with the 5Mb dinucs)
data <- read.delim(list.files(path, "genome_dinucleotide.txt", full.names = TRUE))
normal <- read.delim(list.files(healthy_path, "genome_dinucleotide.txt", full.names = TRUE))
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data$sample, ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
### Sum frequencies across A/T and G/C combinations (LFS)
data$context <- ifelse(data$context %in% c("AA", "AT", "TA", "TT"), "at",
                       ifelse(data$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
data <- data[complete.cases(data), ]
data_contexts <- aggregate(.~context+sample, data, sum)

data_at <- data_contexts[data_contexts$context == "at",]
row.names(data_at) <- data_at$sample
colnames(data_at) <- paste0(colnames(data_at), "_AT")
data_at <- data_at[, colnames(data_at) %like% "X"]

data_gc <- data_contexts[data_contexts$context == "gc",]
row.names(data_gc) <- data_gc$sample
colnames(data_gc) <- paste0(colnames(data_gc), "_GC")
data_gc <- data_gc[, colnames(data_gc) %like% "X"]

data <- bind_cols(data_at, data_gc)
data <- data[data_samples$sWGS, ]

### (Normal)
normal$context <- ifelse(normal$context %in% c("AA", "AT", "TA", "TT"), "at",
                              ifelse(normal$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
normal <- normal[complete.cases(normal), ]
normal_contexts <- aggregate(.~context+sample, normal, sum)

normal_at <- normal_contexts[normal_contexts$context == "at",]
row.names(normal_at) <- normal_at$sample
colnames(normal_at) <- paste0(colnames(normal_at), "_AT")
normal_at <- normal_at[, colnames(normal_at) %like% "X"]

normal_gc <- normal_contexts[normal_contexts$context == "gc",]
row.names(normal_gc) <- normal_gc$sample
colnames(normal_gc) <- paste0(colnames(normal_gc), "_GC")
normal_gc <- normal_gc[, colnames(normal_gc) %like% "X"]

normal <- bind_cols(normal_at, normal_gc)

### Run classifiers
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/runner.R")
