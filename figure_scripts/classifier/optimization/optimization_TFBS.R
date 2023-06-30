library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
path <- "data/griffin/TFBS"
healthy_path <- "hbc/griffin/TFBS"
outdir <- ""
title <- "TFBS"
title2 <- "Nuc (TFBS)"

### Import data 
data <- list.files(path, "features", full.names = TRUE)
normal <- list.files(healthy_path, "features", full.names = TRUE)
data_samples <- read.delim("sample_list.txt")

datalist <- lapply(data, function(x){read.delim(file = x)})
data <- do.call(rbind, datalist)
datalist <- lapply(normal, function(x){read.delim(file = x)})
normal <- do.call(rbind, datalist)

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
row.names(data) <- data$features
data <- data[, data_samples$sWGS]
data <- as.data.frame(t(data))

row.names(normal) <- normal$features
normal <- normal[, -1]
normal <- as.data.frame(t(normal))

### Run classifiers
source("figure_scripts/classifier/optimization/runner.R")
