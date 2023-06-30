library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/end_motifs"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/end_motifs"
butler_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/HBC_Butler/end_motifs"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Import data 
data_end <- read.delim(list.files(path, "CHARM_LFS_genome_end_motifs.txt", full.names = TRUE))
normal_end <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_end_motifs.txt", full.names = TRUE))
butler_end <- read.delim(list.files(butler_path, "HBC_Butler_end_motifs.txt", full.names = TRUE))
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_end), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
data_end[is.na(data_end)] <- 0
row.names(data_end) <- data_end$motif
data_end <- data_end[, colnames(data_end) %in% data_samples$sWGS]
data_end <- as.data.frame(t(data_end))
data_end <- data_end[data_samples$sWGS, ]

normal_end[is.na(normal_end)] <- 0
row.names(normal_end) <- normal_end$motif
normal_end <- normal_end[, -1]
normal_end <- as.data.frame(t(normal_end))

butler_end[is.na(butler_end)] <- 0
row.names(butler_end) <- butler_end$motif
butler_end <- butler_end[, -1]
butler_end <- as.data.frame(t(butler_end))

### Make Healthy vs LFS Negative
Data <- data_end[row.names(data_end) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_end)
Data$sample <- row.names(Data)

Data_butler <- butler_end[complete.cases(butler_end), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_neg <- All.kFold
outcomes_neg_but <- But.kfold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_end_motifs_status.rds"))
saveRDS(outcomes_neg_but, file.path(outdir, "outcomes_end_motifs_status_butler.rds"))

### Make LFS positive vs negative
Data <- data_end[row.names(data_end) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_lfs <- All.kFold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_end_motifs_lfs.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_end_motifs_lfs_butler.rds"))
