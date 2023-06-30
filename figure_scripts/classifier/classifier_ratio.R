library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
butler <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/HBC_Butler"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Import data (Starting with the 5Mb ratios)
data_ratio <- read.delim(file.path(path, "fragmentomics", "CHARM_LFS_ratio_5Mb.txt"))
normal_ratio <- read.delim(file.path(healthy, "fragmentomics", "TGL49_HBC_ratio_5Mb.txt"))
butler_ratio <- read.delim(file.path(butler, "fragmentomics", "HBC_Butler_ratio_5Mb.txt"))
data_samples <- read.delim(file.path(path, "samples/sample_list.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_ratio), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
row.names(data_ratio) <- paste0(data_ratio$seqnames, "_", data_ratio$start)
data_ratio <- data_ratio[, data_samples$sWGS]
data_ratio <- as.data.frame(t(data_ratio))
data_ratio <- data_ratio[data_samples$sWGS, ]

row.names(normal_ratio) <- paste0(normal_ratio$seqnames, "_", normal_ratio$start)
normal_ratio <- normal_ratio[, !(colnames(normal_ratio) %in% c("seqnames", "arm", "start", "end"))]
normal_ratio <- as.data.frame(t(normal_ratio))

row.names(butler_ratio) <- paste0(butler_ratio$seqnames, "_", butler_ratio$start)
butler_ratio <- butler_ratio[, !(colnames(butler_ratio) %in% c("seqnames", "arm", "start", "end"))]
butler_ratio <- as.data.frame(t(butler_ratio))

### Make Healthy vs LFS Negative
Data <- data_ratio[row.names(data_ratio) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_ratio)
Data$sample <- row.names(Data)

Data_butler <- butler_ratio[complete.cases(butler_ratio), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "svmRadial"
source(class)
outcomes_neg <- All.kFold
outcomes_neg_but <- But.kfold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_ratios_status.rds"))
saveRDS(outcomes_neg_but, file.path(outdir, "outcomes_ratio_status_butler.rds"))

### Make LFS positive vs negative
Data <- data_ratio[row.names(data_ratio) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "svmRadial"
source(class)
outcomes_lfs <- All.kFold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_ratios_lfs.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_ratio_lfs_butler.rds"))
