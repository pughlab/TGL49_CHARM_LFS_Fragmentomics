library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization"
title <- "size"
title2 <- "Fragment Length"

### Import data (Starting with the 5Mb sizes)
data_samples <- read.delim(file.path(path, "samples/sample_list.txt"))
data <- read.delim(file.path(path, "insert_size", "CHARM_LFS_fragment_freq.txt"))
normal <- read.delim(file.path(healthy, "insert_size", "TGL49_HBC_fragment_freq.txt"))

data_prop <- read.delim(file.path(path, "insert_size", "CHARM_LFS_proportions.txt"))
normal_prop <- read.delim(file.path(healthy, "insert_size", "TGL49_HBC_proportions.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format frequency table
row.names(data) <- data$length
data <- data[, data_samples$sWGS]
data <- as.data.frame(t(data))

row.names(normal) <- normal$length
normal <- normal[, -1]
normal <- as.data.frame(t(normal))

### Append fragment proportions to frequency tables
data_prop <- data_prop[data_prop$sample %in% data_samples$sWGS, ]
data_prop <- data_prop[order(factor(data_prop$sample, levels = row.names(data))), ]
data$prop <- data_prop$P.20_150.
data <- data[data_samples$sWGS, ]

normal_prop <- normal_prop[normal_prop$sample %in% row.names(normal), ]
normal_prop <- normal_prop[order(factor(normal_prop$sample, levels = row.names(normal))), ]
normal$prop <- normal_prop$P.20_150.

### Run classifiers
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/runner.R")
