library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(data.table)

### Set variables
source("figure_scripts/classifier_split/FuncClassifier.R")
class <- "figure_scripts/classifier_split/classifier.R"
path <- "data/classifier_split"
outdir <- ""

### Read in sample data
data_samples <- read.delim(file.path(outdir, "samples_split.txt"))

### Remove failed data_samples
samples_hbc_test <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "test", ]
samples_hbc_val <- data_samples[data_samples$cancer_status == "healthy" & data_samples$split == "validation", ]
samples_neg_test <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "test", ]
samples_neg_val <- data_samples[data_samples$cancer_status == "negative" & data_samples$split == "validation", ]
samples_pos_test <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "test", ]
samples_pos_val <- data_samples[data_samples$cancer_status == "positive" & data_samples$split == "validation", ]

### List classifiers
class_list <- list.files(path, "lfs.rds", full.names = TRUE)
class_list <- class_list[!grepl("integrated", class_list)]

names <- list.files(path, "lfs.rds", full.names = FALSE)
names <- names[!grepl("integrated", names)]
names <- gsub("outcomes_", "", names)
names <- gsub("_lfs.rds", "", names)

file <- class_list[[1]]
name <- names[[1]]
data <- readRDS(file)
data <- unlist(data, recursive = FALSE)
data <- bind_rows(data)
data$sample <- make.names(data$sample, unique=TRUE)
data <- data[, c("positive", "ActualClass", "sample")]
colnames(data) <- c(name, "ActualClass", "sample")

### Read in and merge data
for (i in c(2:length(class_list))) {
  ### Set variables
  file <- class_list[[i]]
  name <- names[[i]]
  
  ### Read in data and format
  a <- readRDS(file)
  a <- unlist(a, recursive = FALSE)
  a <- bind_rows(a)
  a$sample <- make.names(a$sample, unique=TRUE)
  a <- a[, c("positive", "ActualClass", "sample")]
  colnames(a) <- c(name, "ActualClass", "sample")
  
  ### Merge data
  data <- merge(data, a, by = c("ActualClass", "sample"), all = TRUE)
}

### Format data by aggregating by sample
data <- data[complete.cases(data), ]
data$sample <- gsub("(_WG).*", "\\1", data$sample)

Data <- aggregate(data[, 3:ncol(data)], list(data$sample, data$ActualClass), mean)
row.names(Data) <- Data$Group.1
Data <- Data[, names]
Data_validation <- Data[c(samples_neg_val$sWGS, samples_pos_val$sWGS), ]
Data_validation <- Data_validation[complete.cases(Data_validation), ]

Data <- Data[c(samples_neg_test$sWGS, samples_pos_test$sWGS), ]
Data$sample <- row.names(Data)
Data <- Data[complete.cases(Data), ]

y <- c(rep("negative", sum(row.names(Data) %in% samples_neg_test$sWGS)))
y <- c(y, rep("positive", sum(!(row.names(Data) %in% samples_neg_test$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

### Run classifier
algorithm <- "rf"
source(class)
outcomes_neg_test <- All.kFold
outcomes_neg_val <- Val.kfold
saveRDS(outcomes_neg_test, file.path(outdir, "outcomes_integrated_lfs_test.rds"))
saveRDS(outcomes_neg_val, file.path(outdir, "outcomes_integrated_lfs_val.rds"))

