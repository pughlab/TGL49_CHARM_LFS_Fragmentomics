library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split"

### Import data (Starting with the 5Mb ratios)
data_ratio <- read.delim(file.path(path, "fragmentomics", "CHARM_LFS_ratio_5Mb.txt"))
normal_ratio <- read.delim(file.path(healthy, "fragmentomics", "TGL49_HBC_ratio_5Mb.txt"))
data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### Format table for classifier
row.names(data_ratio) <- paste0(data_ratio$seqnames, "_", data_ratio$start)
data_ratio <- data_ratio[, colnames(data_ratio) %in% data_samples$sWGS]
data_ratio <- as.data.frame(t(data_ratio))

row.names(normal_ratio) <- paste0(normal_ratio$seqnames, "_", normal_ratio$start)
normal_ratio <- normal_ratio[, !(colnames(normal_ratio) %in% c("seqnames", "arm", "start", "end"))]
normal_ratio <- as.data.frame(t(normal_ratio))

### Make Healthy vs LFS Negative
Data <- bind_rows(data_ratio[samples_neg_test$sWGS, ], 
                  normal_ratio[samples_hbc_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_ratio[samples_neg_val$sWGS, ], 
                             normal_ratio[samples_hbc_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))


algorithm <- "svmRadial"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_ratio_status_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_ratio_status_val.rds"))

### Make LFS positive vs negative
Data <- bind_rows(data_ratio[samples_neg_test$sWGS, ], 
                  data_ratio[samples_pos_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_ratio[samples_neg_val$sWGS, ], 
                             data_ratio[samples_pos_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "svmRadial"
source(class)
outcomes_lfs_test <- All.kFold
outcomes_lfs_val <- Val.kfold
saveRDS(outcomes_lfs_test, file.path(outdir, "outcomes_ratio_lfs_test.rds"))
saveRDS(outcomes_lfs_val, file.path(outdir, "outcomes_ratio_lfs_val.rds"))
