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

### Import data (Starting with the 5Mb sizes)
data_freq <- read.delim(file.path(path, "insert_size", "CHARM_LFS_fragment_freq.txt"))
normal_freq <- read.delim(file.path(healthy, "insert_size", "TGL49_HBC_fragment_freq.txt"))

data_prop <- read.delim(file.path(path, "insert_size", "CHARM_LFS_proportions.txt"))
normal_prop <- read.delim(file.path(healthy, "insert_size", "TGL49_HBC_proportions.txt"))

data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### Format frequency table
row.names(data_freq) <- data_freq$length
data_freq <- data_freq[, colnames(data_freq) %in% data_samples$sWGS]
data_freq <- as.data.frame(t(data_freq))

row.names(normal_freq) <- normal_freq$length
normal_freq <- normal_freq[, -1]
normal_freq <- as.data.frame(t(normal_freq))

### Append fragment proportions to frequency tables
data_prop <- data_prop[data_prop$sample %in% data_samples$sWGS, ]
data_prop <- data_prop[order(factor(data_prop$sample, levels = row.names(data_freq))), ]
data_freq$prop <- data_prop$P.20_150.

normal_prop <- normal_prop[normal_prop$sample %in% row.names(normal_freq), ]
normal_prop <- normal_prop[order(factor(normal_prop$sample, levels = row.names(normal_freq))), ]
normal_freq$prop <- normal_prop$P.20_150.

### Make Healthy vs LFS Negative
Data <- bind_rows(data_freq[samples_neg_test$sWGS, ], 
                  normal_freq[samples_hbc_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_freq[samples_neg_val$sWGS, ], 
                             normal_freq[samples_hbc_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "svmRadial"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_insert_size_status_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_insert_size_status_val.rds"))

### Make LFS positive vs negative
Data <- bind_rows(data_freq[samples_neg_test$sWGS, ], 
                  data_freq[samples_pos_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_freq[samples_neg_val$sWGS, ], 
                             data_freq[samples_pos_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "knn"
source(class)
outcomes_lfs_test <- All.kFold
outcomes_lfs_val <- Val.kfold
saveRDS(outcomes_lfs_test, file.path(outdir, "outcomes_insert_size_lfs_test.rds"))
saveRDS(outcomes_lfs_val, file.path(outdir, "outcomes_insert_size_lfs_val.rds"))
