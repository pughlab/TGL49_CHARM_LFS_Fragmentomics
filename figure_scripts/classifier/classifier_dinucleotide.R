library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(data.table)

### Set variables
source("figure_scripts/classifier/FuncClassifier.R")
class <- "figure_scripts/classifier/classifier.R"
path <- "data/dinucleotide"
healthy_path <- "hbc/dinucleotide"
butler_path <- "butler/dinucleotide"
outdir <- ""

### Import data (Starting with the 5Mb dinucs)
data_dinuc <- read.delim(list.files(path, "genome_dinucleotide.txt", full.names = TRUE))
normal_dinuc <- read.delim(list.files(healthy_path, "genome_dinucleotide.txt", full.names = TRUE))
butler_dinuc <- read.delim(list.files(butler_path, "genome_dinucleotide.txt", full.names = TRUE))
data_samples <- read.delim("sample_list.txt")

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

### (Butler)
butler_dinuc$context <- ifelse(butler_dinuc$context %in% c("AA", "AT", "TA", "TT"), "at",
                               ifelse(butler_dinuc$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
butler_dinuc <- butler_dinuc[complete.cases(butler_dinuc), ]
butler_dinuc_contexts <- aggregate(.~context+sample, butler_dinuc, sum)

butler_dinuc_at <- butler_dinuc_contexts[butler_dinuc_contexts$context == "at",]
row.names(butler_dinuc_at) <- butler_dinuc_at$sample
colnames(butler_dinuc_at) <- paste0(colnames(butler_dinuc_at), "_AT")
butler_dinuc_at <- butler_dinuc_at[, colnames(butler_dinuc_at) %like% "X"]

butler_dinuc_gc <- butler_dinuc_contexts[butler_dinuc_contexts$context == "gc",]
row.names(butler_dinuc_gc) <- butler_dinuc_gc$sample
colnames(butler_dinuc_gc) <- paste0(colnames(butler_dinuc_gc), "_GC")
butler_dinuc_gc <- butler_dinuc_gc[, colnames(butler_dinuc_gc) %like% "X"]

butler_dinuc_contexts <- bind_cols(butler_dinuc_at, butler_dinuc_gc)

### Make Healthy vs LFS Negative
Data <- data_contexts[row.names(data_contexts) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_dinuc_contexts)
Data$sample <- row.names(Data)

Data_butler <- butler_dinuc_contexts[complete.cases(butler_dinuc_contexts), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_neg <- All.kFold
outcomes_neg_but <- But.kfold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_dinucleotide_status.rds"))
saveRDS(outcomes_neg_but, file.path(outdir, "outcomes_dinucleotide_status_butler.rds"))

### Make LFS positive vs negative
Data <- data_contexts[row.names(data_contexts) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "glm"
source(class)
outcomes_lfs <- All.kFold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_dinucleotide_lfs.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_dinucleotide_lfs_butler.rds"))

