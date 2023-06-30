library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(multiROC)

### Set variables
path <- "data/griffin"
healthy_path <- "hbc/griffin"
outdir <- ""
title <- "TCGA"
title2 <- "Nuc (TCGA/DHS)"

sites <- c("DHS", "TCGA")
data <- c()
normal <- c()

for (site in sites) {
  files_griffin <- list.files(file.path(path, site), "features", full.names = TRUE)
  files_normal <- list.files(file.path(healthy_path, site), "features", full.names = TRUE)
  
  datalist <- lapply(files_griffin, function(x){read.delim(file = x)})
  datalist <- do.call(rbind, datalist)
  data <- rbind(data, datalist)
  
  datalist <- lapply(files_normal, function(x){read.delim(file = x)})
  datalist <- do.call(rbind, datalist)
  normal <- rbind(normal, datalist)
}
### Import data 
data_samples <- read.delim("sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
row.names(data) <- data$features
data <- data[, colnames(data) %in% data_samples$sWGS]
data <- data[, data_samples$sWGS]
data <- as.data.frame(t(data))

row.names(normal) <- normal$features
normal <- normal[, -1]
normal <- as.data.frame(t(normal))

### Run classifiers
source("figures_scripts/classifier/optimization/runner.R")
