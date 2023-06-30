library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/nucleosome_peaks"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/nucleosome_peaks"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split"

### Import data 
data_peak <- read.delim(list.files(path, "CHARM_LFS_genome_peak_distances.txt", full.names = TRUE))
normal_peak <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_peak_distances.txt", full.names = TRUE))
data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### Format table for classifier
data_peak <- data_peak[data_peak$sample %in% data_samples$sWGS, !(colnames(data_peak) %in% c("chr"))]
names <- unique(data_peak$sample)
data_peak_sum <- aggregate(. ~ sample, data_peak, sum)
row.names(data_peak_sum) <- data_peak_sum$sample
data_peak_sum <- data_peak_sum[,-1]
colnames(data_peak_sum) <- c(-1000:1000)
data_peak_sum <- data_peak_sum[, colnames(data_peak_sum) %in% c(-300:300)]  
data_peak_sum <- data_peak_sum/rowSums(data_peak_sum)

normal_peak <- normal_peak[, !(colnames(normal_peak) %in% c("chr"))]
names <- unique(normal_peak$sample)
normal_peak_sum <- aggregate(. ~ sample, normal_peak, sum)
row.names(normal_peak_sum) <- normal_peak_sum$sample
normal_peak_sum <- normal_peak_sum[,-1]
colnames(normal_peak_sum) <- c(-1000:1000)
normal_peak_sum <- normal_peak_sum[, colnames(normal_peak_sum) %in% c(-300:300)]  
normal_peak_sum <- normal_peak_sum/rowSums(normal_peak_sum)

### Make Healthy vs LFS Negative
Data <- bind_rows(data_peak_sum[samples_neg_test$sWGS, ], 
                  normal_peak_sum[samples_hbc_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_peak_sum[samples_neg_val$sWGS, ], 
                             normal_peak_sum[samples_hbc_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_peaks_status_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_peaks_status_val.rds"))

### Make LFS positive vs negative
Data <- bind_rows(data_peak_sum[samples_neg_test$sWGS, ], 
                  data_peak_sum[samples_pos_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_peak_sum[samples_neg_val$sWGS, ], 
                             data_peak_sum[samples_pos_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "knn"
source(class)
outcomes_lfs_test <- All.kFold
outcomes_lfs_val <- Val.kfold
saveRDS(outcomes_lfs_test, file.path(outdir, "outcomes_peaks_lfs_test.rds"))
saveRDS(outcomes_lfs_val, file.path(outdir, "outcomes_peaks_lfs_val.rds"))
