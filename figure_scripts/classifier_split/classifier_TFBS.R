library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/griffin_all/TFBS"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/griffin_all/TFBS"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split"

### Import data 
data_griffin <- list.files(path, "features", full.names = TRUE)
normal_griffin <- list.files(healthy_path, "features", full.names = TRUE)

datalist <- lapply(data_griffin, function(x){read.delim(file = x)})
data_griffin <- do.call(rbind, datalist)
datalist <- lapply(normal_griffin, function(x){read.delim(file = x)})
normal_griffin <- do.call(rbind, datalist)

data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### Format table for classifier
row.names(data_griffin) <- data_griffin$features
data_griffin <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS]
data_griffin <- as.data.frame(t(data_griffin))

row.names(normal_griffin) <- normal_griffin$features
normal_griffin <- normal_griffin[, -1]
normal_griffin <- as.data.frame(t(normal_griffin))

### Make Healthy vs LFS Negative
Data <- bind_rows(data_griffin[samples_neg_test$sWGS, ], 
                  normal_griffin[samples_hbc_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_griffin[samples_neg_val$sWGS, ], 
                             normal_griffin[samples_hbc_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]


y <- c(rep("positive", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_TFBS_status_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_TFBS_status_val.rds"))

### Make LFS positive vs negative
Data <- bind_rows(data_griffin[samples_neg_test$sWGS, ], 
                  data_griffin[samples_pos_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_griffin[samples_neg_val$sWGS, ], 
                             data_griffin[samples_pos_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_lfs_test <- All.kFold
outcomes_lfs_val <- Val.kfold
saveRDS(outcomes_lfs_test, file.path(outdir, "outcomes_TFBS_lfs_test.rds"))
saveRDS(outcomes_lfs_val, file.path(outdir, "outcomes_TFBS_lfs_val.rds"))
