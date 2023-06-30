library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/breakpoint"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/breakpoint"
butler_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/HBC_Butler/breakpoint"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Import data 
data <- read.delim(list.files(path, "CHARM_LFS_genome_breakpoint_ratio.txt", full.names = TRUE))
normal <- read.delim(list.files(healthy_path, "CHARM_HBC_genome_breakpoint_ratio.txt", full.names = TRUE))
butler <- read.delim(list.files(butler_path, "HBC_Butler_output_breakpoint_ratio.txt", full.names = TRUE))
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% data$sample, ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
data[is.na(data)] <- 0
data <- reshape2::melt(data, id = c("nucleotide", "sample"))
data$context <- paste0(data$nucleotide, "_", data$variable)
data <- reshape2::dcast(data[, c("sample", "value", "context")], context ~ sample)

row.names(data) <- data$context
data <- data[, colnames(data) %in% data_samples$sWGS]
data <- as.data.frame(t(data))
data <- data[data_samples$sWGS, ]

normal[is.na(normal)] <- 0
normal <- reshape2::melt(normal, id = c("nucleotide", "sample"))
normal$context <- paste0(normal$nucleotide, "_", normal$variable)
normal <- reshape2::dcast(normal[, c("sample", "value", "context")], context ~ sample)

row.names(normal) <- normal$context
normal <- normal[, -1]
normal <- as.data.frame(t(normal))

butler[is.na(butler)] <- 0
butler <- reshape2::melt(butler, id = c("nucleotide", "sample"))
butler$context <- paste0(butler$nucleotide, "_", butler$variable)
butler <- reshape2::dcast(butler[, c("sample", "value", "context")], context ~ sample)

row.names(butler) <- butler$context
butler <- butler[, -1]
butler <- as.data.frame(t(butler))

### Make Healthy vs LFS Negative
Data <- data[row.names(data) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal)
Data$sample <- row.names(Data)

Data_butler <- butler[complete.cases(butler), ]

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "svmRadial"
source(class)
outcomes_neg <- All.kFold
outcomes_neg_but <- But.kfold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_breakpoint_status.rds"))
saveRDS(outcomes_neg_but, file.path(outdir, "outcomes_breakpoint_status_butler.rds"))

### Make LFS positive vs negative
Data <- data[row.names(data) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_lfs <- All.kFold
outcomes_lfs_but <- But.kfold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_breakpoint_lfs.rds"))
saveRDS(outcomes_lfs_but, file.path(outdir, "outcomes_breakpoint_lfs_butler.rds"))
