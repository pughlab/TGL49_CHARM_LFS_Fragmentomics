library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/nucleosome_peaks"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/nucleosome_peaks"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Import data 
data_peak <- read.delim(list.files(path, "CHARM_LFS_genome_peak_distances.txt", full.names = TRUE))
normal_peak <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_peak_distances.txt", full.names = TRUE))
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

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

### Make Healthy vs LFS Negative
Data <- data_peak_sum[row.names(data_peak_sum) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_peak_sum)
Data$sample <- row.names(Data)

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_neg <- All.kFold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_peaks_status.rds"))

### Make LFS positive vs negative
Data <- data_peak_sum[row.names(data_peak_sum) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "knn"
source(class)
outcomes_lfs <- All.kFold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_peaks_lfs.rds"))
