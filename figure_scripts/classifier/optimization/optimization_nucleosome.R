library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
path <- "data/nucleosome_peaks"
healthy_path <- "hbc/nucleosome_peaks"
outdir <- ""
title <- "nucleosome"
title2 <- "Nuc Dist"

### Import data 
data <- read.delim(list.files(path, "CHARM_LFS_genome_peak_distances.txt", full.names = TRUE))
normal <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_peak_distances.txt", full.names = TRUE))
data_samples <- read.delim("sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data$sample, ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
data <- data[data$sample %in% data_samples$sWGS, !(colnames(data) %in% c("chr"))]
names <- unique(data$sample)
data_sum <- aggregate(. ~ sample, data, sum)
row.names(data_sum) <- data_sum$sample
data_sum <- data_sum[,-1]
colnames(data_sum) <- c(-1000:1000)
data_sum <- data_sum[, colnames(data_sum) %in% c(-300:300)]  
data_sum <- data_sum/rowSums(data_sum)
data <- data_sum[data_samples$sWGS,]

normal <- normal[, !(colnames(normal) %in% c("chr"))]
names <- unique(normal$sample)
normal_sum <- aggregate(. ~ sample, normal, sum)
row.names(normal_sum) <- normal_sum$sample
normal_sum <- normal_sum[,-1]
colnames(normal_sum) <- c(-1000:1000)
normal_sum <- normal_sum[, colnames(normal_sum) %in% c(-300:300)]  
normal_sum <- normal_sum/rowSums(normal_sum)
normal <- normal_sum

### Run classifiers
source("figure_scripts/classifier/optimization/runner.R")
