library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/griffin_all/TFBS"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/griffin_all/TFBS"
butler_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/HBC_Butler/griffin_all/TFBS"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Import data 
data_griffin <- list.files(path, "features", full.names = TRUE)
normal_griffin <- list.files(healthy_path, "features", full.names = TRUE)
butler_griffin <- list.files(butler_path, "features", full.names = TRUE)

data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

datalist <- lapply(data_griffin, function(x){read.delim(file = x)})
data_griffin <- do.call(rbind, datalist)
datalist <- lapply(normal_griffin, function(x){read.delim(file = x)})
normal_griffin <- do.call(rbind, datalist)
datalist <- lapply(butler_griffin, function(x){read.delim(file = x)})
butler_griffin <- do.call(rbind, datalist)

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_griffin), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
row.names(data_griffin) <- data_griffin$features
data_griffin <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS]
data_griffin <- as.data.frame(t(data_griffin))
data_griffin <- data_griffin[data_samples$sWGS, ]

row.names(normal_griffin) <- normal_griffin$features
normal_griffin <- normal_griffin[, -1]
normal_griffin <- as.data.frame(t(normal_griffin))

row.names(butler_griffin) <- butler_griffin$features
butler_griffin <- butler_griffin[, -1]
butler_griffin <- as.data.frame(t(butler_griffin))

### Make Healthy vs LFS Negative
Data <- data_griffin[row.names(data_griffin) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_griffin)
Data$sample <- row.names(Data)

Data_butler <- butler_griffin[complete.cases(butler_griffin), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_neg <- All.kFold
outcomes_neg_but <- But.kfold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_TFBS_status.rds"))
saveRDS(outcomes_neg_but, file.path(outdir, "outcomes_TFBS_status_butler.rds"))

### Make LFS positive vs negative
Data <- data_griffin[row.names(data_griffin) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_lfs <- All.kFold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_TFBS_lfs.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_TFBS_lfs_butler.rds"))
