library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### List classifiers
class_list <- list.files(path, "status.rds", full.names = TRUE)
class_list <- class_list[!grepl("integrated", class_list)]
class_list <- class_list[!grepl("end_motif", class_list)]

names <- list.files(path, "status.rds", full.names = FALSE)
names <- names[!grepl("integrated", names)]
names <- names[!grepl("end_motif", names)]
names <- gsub("outcomes_", "", names)
names <- gsub("_status.rds", "", names)

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
Data$sample <- Data$Group.1
y <- Data$Group.2
y <- factor(y, levels = c("negative", "positive"))
Data <- Data[, !(colnames(Data) %in% c("Group.1", "Group.2"))]

### Run classifier
algorithm <- "rf"
source(class)
outcomes_lfs <- All.kFold
saveRDS(outcomes_lfs, file.path(path, "outcomes_integrated_status.rds"))
