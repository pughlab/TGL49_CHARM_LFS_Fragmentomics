library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("figure_scripts/classifier/FuncClassifier.R")
class <- "figure_scripts/classifier/classifier.R"
path <- "data/nucleosome_peaks"
healthy_path <- "hbc/nucleosome_peaks"
butler_path <- "butler/nucleosome_peaks"
outdir <- ""

### Import data 
data_peak <- read.delim(list.files(path, "CHARM_LFS_genome_peak_distances.txt", full.names = TRUE))
normal_peak <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_peak_distances.txt", full.names = TRUE))
butler_peak <- read.delim(list.files(butler_path, "HBC_Butler_genome_peak_distances.txt", full.names = TRUE))
data_samples <- read.delim("sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data_peak$sample, ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
data_peak <- data_peak[data_peak$sample %in% data_samples$sWGS, !(colnames(data_peak) %in% c("chr"))]
names <- unique(data_peak$sample)
data_peak_sum <- aggregate(. ~ sample, data_peak, sum)
row.names(data_peak_sum) <- data_peak_sum$sample
data_peak_sum <- data_peak_sum[,-1]
colnames(data_peak_sum) <- c(-1000:1000)
data_peak_sum <- data_peak_sum[, colnames(data_peak_sum) %in% c(-300:300)]  
data_peak_sum <- data_peak_sum/rowSums(data_peak_sum)
data_peak_sum <- data_peak_sum[data_samples$sWGS,]

normal_peak <- normal_peak[, !(colnames(normal_peak) %in% c("chr"))]
names <- unique(normal_peak$sample)
normal_peak_sum <- aggregate(. ~ sample, normal_peak, sum)
row.names(normal_peak_sum) <- normal_peak_sum$sample
normal_peak_sum <- normal_peak_sum[,-1]
colnames(normal_peak_sum) <- c(-1000:1000)
normal_peak_sum <- normal_peak_sum[, colnames(normal_peak_sum) %in% c(-300:300)]  
normal_peak_sum <- normal_peak_sum/rowSums(normal_peak_sum)

butler_peak <- butler_peak[, !(colnames(butler_peak) %in% c("chr"))]
names <- unique(butler_peak$sample)
butler_peak_sum <- aggregate(. ~ sample, butler_peak, sum)
row.names(butler_peak_sum) <- butler_peak_sum$sample
butler_peak_sum <- butler_peak_sum[,-1]
colnames(butler_peak_sum) <- c(-1000:1000)
butler_peak_sum <- butler_peak_sum[, colnames(butler_peak_sum) %in% c(-300:300)]  
butler_peak_sum <- butler_peak_sum/rowSums(butler_peak_sum)

### Make Healthy vs LFS Negative
Data <- data_peak_sum[row.names(data_peak_sum) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_peak_sum)
Data$sample <- row.names(Data)

Data_butler <- butler_peak_sum[complete.cases(butler_peak_sum), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_neg <- All.kFold
outcomes_neg_but <- But.kfold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_peaks_status.rds"))
saveRDS(outcomes_neg_but, file.path(outdir, "outcomes_peaks_status_butler.rds"))

### Make LFS positive vs negative
Data <- data_peak_sum[row.names(data_peak_sum) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "knn"
source(class)
outcomes_lfs <- All.kFold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_peaks_lfs.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_peaks_lfs_butler.rds"))
