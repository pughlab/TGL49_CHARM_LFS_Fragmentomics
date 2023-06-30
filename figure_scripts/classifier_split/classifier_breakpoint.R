library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/FuncClassifier.R")
class <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/breakpoint"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/breakpoint"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split"

### Import data 
data <- read.delim(list.files(path, "CHARM_LFS_genome_breakpoint_ratio.txt", full.names = TRUE))
normal <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_breakpoint_ratio.txt", full.names = TRUE))
data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### Format table for classifier
data[is.na(data)] <- 0
data <- reshape2::melt(data, id = c("nucleotide", "sample"))
data$context <- paste0(data$nucleotide, "_", data$variable)
data <- reshape2::dcast(data[, c("sample", "value", "context")], context ~ sample)

row.names(data) <- data$context
data <- data[, colnames(data) %in% data_samples$sWGS]
data <- as.data.frame(t(data))

normal[is.na(normal)] <- 0
normal <- reshape2::melt(normal, id = c("nucleotide", "sample"))
normal$context <- paste0(normal$nucleotide, "_", normal$variable)
normal <- reshape2::dcast(normal[, c("sample", "value", "context")], context ~ sample)

row.names(normal) <- normal$context
normal <- normal[, -1]
normal <- as.data.frame(t(normal))

### Make Healthy vs LFS Negative
Data <- bind_rows(data[samples_neg_test$sWGS, ], 
                  normal[samples_hbc_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data[samples_neg_val$sWGS, ], 
                             normal[samples_hbc_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "svmRadial"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_breakpoint_status_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_breakpoint_status_val.rds"))

### Make LFS positive vs negative
Data <- bind_rows(data[samples_neg_test$sWGS, ], 
                  data[samples_pos_test$sWGS, ])
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

Data_validation <- bind_rows(data[samples_neg_val$sWGS, ], 
                             data[samples_pos_val$sWGS, ])
Data_validation <- Data_validation[complete.cases(Data_validation), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_lfs_test <- All.kFold
outcomes_lfs_val <- Val.kfold
saveRDS(outcomes_lfs_test, file.path(outdir, "outcomes_breakpoint_lfs_test.rds"))
saveRDS(outcomes_lfs_val, file.path(outdir, "outcomes_breakpoint_lfs_val.rds"))
