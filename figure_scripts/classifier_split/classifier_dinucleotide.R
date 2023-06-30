library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(data.table)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/dinucleotide"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/dinucleotide"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split"

### Import data (Starting with the 5Mb dinucs)
data_dinuc <- read.delim(list.files(path, "genome_dinucleotide.txt", full.names = TRUE))
normal_dinuc <- read.delim(list.files(healthy_path, "genome_dinucleotide.txt", full.names = TRUE))
data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### Format table for classifier
### Sum frequencies across A/T and G/C combinations (LFS)
data_dinuc <- data_dinuc[data_dinuc$sample %in% data_samples$sWGS, ]
data_dinuc$context <- ifelse(data_dinuc$context %in% c("AA", "AT", "TA", "TT"), "at",
                       ifelse(data_dinuc$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
data_dinuc <- data_dinuc[complete.cases(data_dinuc), ]
data_contexts <- aggregate(.~context+sample, data_dinuc, sum)

data_at <- data_contexts[data_contexts$context == "at",]
row.names(data_at) <- data_at$sample
colnames(data_at) <- paste0(colnames(data_at), "_AT")
data_at <- data_at[, colnames(data_at) %like% "X"]

data_gc <- data_contexts[data_contexts$context == "gc",]
row.names(data_gc) <- data_gc$sample
colnames(data_gc) <- paste0(colnames(data_gc), "_GC")
data_gc <- data_gc[, colnames(data_gc) %like% "X"]

data_contexts <- bind_cols(data_at, data_gc)

### (Normal)
normal_dinuc$context <- ifelse(normal_dinuc$context %in% c("AA", "AT", "TA", "TT"), "at",
                              ifelse(normal_dinuc$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
normal_dinuc <- normal_dinuc[complete.cases(normal_dinuc), ]
normal_dinuc_contexts <- aggregate(.~context+sample, normal_dinuc, sum)

normal_dinuc_at <- normal_dinuc_contexts[normal_dinuc_contexts$context == "at",]
row.names(normal_dinuc_at) <- normal_dinuc_at$sample
colnames(normal_dinuc_at) <- paste0(colnames(normal_dinuc_at), "_AT")
normal_dinuc_at <- normal_dinuc_at[, colnames(normal_dinuc_at) %like% "X"]

normal_dinuc_gc <- normal_dinuc_contexts[normal_dinuc_contexts$context == "gc",]
row.names(normal_dinuc_gc) <- normal_dinuc_gc$sample
colnames(normal_dinuc_gc) <- paste0(colnames(normal_dinuc_gc), "_GC")
normal_dinuc_gc <- normal_dinuc_gc[, colnames(normal_dinuc_gc) %like% "X"]

normal_dinuc_contexts <- bind_cols(normal_dinuc_at, normal_dinuc_gc)

### Make Healthy vs LFS Negative
Data <- bind_rows(data_contexts[samples_neg_test$sWGS, ], 
                  normal_dinuc_contexts[samples_hbc_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_contexts[samples_neg_val$sWGS, ], 
                             normal_dinuc_contexts[samples_hbc_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_dinucleotide_status_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_dinucleotide_status_val.rds"))

### Make LFS positive vs negative
Data <- bind_rows(data_contexts[samples_neg_test$sWGS, ], 
                  data_contexts[samples_pos_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data_contexts[samples_neg_val$sWGS, ], 
                             data_contexts[samples_pos_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "glm"
source(class)
outcomes_lfs_test <- All.kFold
outcomes_lfs_val <- Val.kfold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs_test, file.path(outdir, "outcomes_dinucleotide_lfs_test.rds"))
saveRDS(outcomes_lfs_val, file.path(outdir, "outcomes_dinucleotide_lfs_val.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_dinucleotide_lfs_butler.rds"))
