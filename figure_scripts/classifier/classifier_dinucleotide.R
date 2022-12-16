library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(data.table)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/dinucleotide"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/dinucleotide"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Import data (Starting with the 5Mb dinucs)
data_dinuc <- read.delim(list.files(path, "genome_dinucleotide.txt", full.names = TRUE))
normal_dinuc <- read.delim(list.files(healthy_path, "genome_dinucleotide.txt", full.names = TRUE))
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data_dinuc$sample, ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

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
data_contexts <- data_contexts[data_samples$sWGS, ]

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
Data <- data_contexts[row.names(data_contexts) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_dinuc_contexts)
Data$sample <- row.names(Data)

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_neg <- All.kFold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_dinucleotide_status.rds"))

### Make LFS positive vs negative
Data <- data_contexts[row.names(data_contexts) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "glm"
source(class)
outcomes_lfs <- All.kFold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_dinucleotide_lfs.rds"))
