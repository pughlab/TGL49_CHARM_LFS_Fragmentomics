library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("figure_scripts/classifier/FuncClassifier.R")
class <- "figure_scripts/classifier/classifier.R"
path <- "data/insert_size"
healthy <- "hbc/insert_size"
butler <- "butler/insert_size"
outdir <- ""

### Import data (Starting with the 5Mb sizes)
data_samples <- read.delim("sample_list.txt"))
data_freq <- read.delim(file.path(path, "CHARM_LFS_fragment_freq.txt"))
normal_freq <- read.delim(file.path(healthy, "TGL49_HBC_fragment_freq.txt"))
butler_freq <- read.delim(file.path(butler, "HBC_Butler_fragment_freq.txt"), check.names = FALSE)

data_prop <- read.delim(file.path(path, "CHARM_LFS_proportions.txt"))
normal_prop <- read.delim(file.path(healthy, "TGL49_HBC_proportions.txt"))
butler_prop <- read.delim(file.path(butler, "HBC_Butler_proportions.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_freq), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format frequency table
row.names(data_freq) <- data_freq$length
data_freq <- data_freq[, data_samples$sWGS]
data_freq <- as.data.frame(t(data_freq))

row.names(normal_freq) <- normal_freq$length
normal_freq <- normal_freq[, -1]
normal_freq <- as.data.frame(t(normal_freq))

row.names(butler_freq) <- butler_freq$length
butler_freq <- butler_freq[, -1]
butler_freq <- as.data.frame(t(butler_freq))

### Append fragment proportions to frequency tables
data_prop <- data_prop[data_prop$sample %in% data_samples$sWGS, ]
data_prop <- data_prop[order(factor(data_prop$sample, levels = row.names(data_freq))), ]
data_freq$prop <- data_prop$P.20_150.
data_freq <- data_freq[data_samples$sWGS, ]

normal_prop <- normal_prop[normal_prop$sample %in% row.names(normal_freq), ]
normal_prop <- normal_prop[order(factor(normal_prop$sample, levels = row.names(normal_freq))), ]
normal_freq$prop <- normal_prop$P.20_150.

butler_prop <- butler_prop[butler_prop$sample %in% row.names(butler_freq), ]
butler_prop <- butler_prop[order(factor(butler_prop$sample, levels = row.names(butler_freq))), ]
butler_freq$prop <- butler_prop$P.20_150.

### Make Healthy vs LFS Negative
Data <- data_freq[row.names(data_freq) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_freq)
Data$sample <- row.names(Data)

Data_butler <- butler_freq[complete.cases(butler_freq), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "svmRadial"
source(class)
outcomes_neg <- All.kFold
outcomes_neg_but <- But.kfold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_insert_size_status.rds"))
saveRDS(outcomes_neg_but, file.path(outdir, "outcomes_insert_size_status_butler.rds"))

### Make LFS positive vs negative
Data <- data_freq[row.names(data_freq) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "knn"
source(class)
outcomes_lfs <- All.kFold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_insert_size_lfs.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_insert_size_lfs_butler.rds"))
