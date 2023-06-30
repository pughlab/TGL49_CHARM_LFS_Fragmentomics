library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/end_motifs"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/end_motifs"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split"

### Import data 
data_end <- read.delim(list.files(path, "CHARM_LFS_genome_end_motifs.txt", full.names = TRUE))
normal_end <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_end_motifs.txt", full.names = TRUE))
data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### Format table for classifier
data_end[is.na(data_end)] <- 0
row.names(data_end) <- data_end$motif
data_end <- data_end[, colnames(data_end) %in% data_samples$sWGS]
data_end <- as.data.frame(t(data_end))

normal_end[is.na(normal_end)] <- 0
row.names(normal_end) <- normal_end$motif
normal_end <- normal_end[, -1]
normal_end <- as.data.frame(t(normal_end))

### Make Healthy vs LFS Negative
Data <- bind_rows(data_end[samples_neg_test$sWGS, ], 
                  normal_end[samples_hbc_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_end[samples_neg_val$sWGS, ], 
                             normal_end[samples_hbc_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_end_motifs_status_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_end_motifs_status_val.rds"))

### Make LFS positive vs negative
Data <- bind_rows(data_end[samples_neg_test$sWGS, ], 
                  data_end[samples_pos_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_end[samples_neg_val$sWGS, ], 
                             data_end[samples_pos_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_lfs_test <- All.kFold
outcomes_lfs_val <- Val.kfold
saveRDS(outcomes_lfs_test, file.path(outdir, "outcomes_end_motifs_lfs_test.rds"))
saveRDS(outcomes_lfs_val, file.path(outdir, "outcomes_end_motifs_lfs_val.rds"))
