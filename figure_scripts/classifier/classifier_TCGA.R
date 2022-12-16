library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/griffin_all"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/griffin_all"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

sites <- c("DHS", "TCGA")
data_griffin <- c()
normal_griffin <- c()

for (site in sites) {
  files_griffin <- list.files(file.path(path, site), "features", full.names = TRUE)
  files_normal <- list.files(file.path(healthy_path, site), "features", full.names = TRUE)
  
  datalist <- lapply(files_griffin, function(x){read.delim(file = x)})
  data <- do.call(rbind, datalist)
  data_griffin <- rbind(data_griffin, data)
  
  datalist <- lapply(files_normal, function(x){read.delim(file = x)})
  data <- do.call(rbind, datalist)
  normal_griffin <- rbind(normal_griffin, data)
}

### Import data 
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_griffin), ]

samples_neg <- data_samples[data_samples$cancer_status == "negative",]
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]

### Format table for classifier
### Sum frequencies across A/T and G/C combinations (LFS)
row.names(data_griffin) <- data_griffin$features
data_griffin <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS]
data_griffin <- as.data.frame(t(data_griffin))
data_griffin <- data_griffin[data_samples$sWGS, ]

row.names(normal_griffin) <- normal_griffin$features
normal_griffin <- normal_griffin[, -1]
normal_griffin <- as.data.frame(t(normal_griffin))

### Make Healthy vs LFS Negative
Data <- data_griffin[row.names(data_griffin) %in% samples_neg$sWGS, ]
Data <- bind_rows(Data, normal_griffin)
Data$sample <- row.names(Data)

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "rf"
source(class)
outcomes_neg <- All.kFold
saveRDS(outcomes_neg, file.path(outdir, "outcomes_TCGA_DHS_status.rds"))

### Make LFS positive vs negative
Data <- data_griffin[row.names(data_griffin) %in% data_samples$sWGS, ]
Data$sample <- row.names(Data)

y <- data_samples$cancer_status
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"
source(class)
outcomes_lfs <- All.kFold
saveRDS(outcomes_lfs, file.path(outdir, "outcomes_TCGA_DHS_lfs.rds"))
