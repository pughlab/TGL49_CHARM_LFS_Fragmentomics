library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(data.table)

### Set variables
source("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/FuncClassifier.R")
class <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/classifier.R"
path <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### List classifiers
class_list <- list.files(path, "status", full.names = TRUE)
class_list <- class_list[grepl(".rds", class_list)]
class_list <- class_list[!grepl("integrated", class_list)]
class_list <- class_list[!grepl("end_motif", class_list)]

names <- list.files(path, "status", full.names = FALSE)
names <- names[grepl(".rds", names)]
names <- names[!grepl("integrated", names)]
names <- names[!grepl("end_motif", names)]
names <- gsub("outcomes_", "", names)
names <- gsub("_status.rds", "", names)
names <- gsub("_status_butler.rds", "", names)

file <- class_list[[1]]
name <- names[[1]]
data <- readRDS(file)
data <- unlist(data, recursive = FALSE)
data <- bind_rows(data)
data$sample <- make.names(data$sample, unique=TRUE)

add_cols <- function(df, cols) {
  add <- cols[!cols %in% names(df)]
  if(length(add) != 0) df[add] <- "negative"
  return(df)
}

data <- add_cols(data, "ActualClass")
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
  a <- add_cols(a, "ActualClass")
  a <- a[, c("positive", "ActualClass", "sample")]
  colnames(a) <- c(name, "ActualClass", "sample")
  
  ### Merge data
  data <- merge(data, a, by = c("ActualClass", "sample"), all = TRUE)
}

### Format data by aggregating by sample
butler <- data[, c("ActualClass", "sample", colnames(data)[colnames(data) %like% ".x"])]
butler <- butler[complete.cases(butler), ]
colnames(butler) <- gsub(".x", "", colnames(butler))
butler$sample <- gsub("(WGS).*", "\\1", butler$sample)

data <- data[, c("ActualClass", "sample", colnames(data)[colnames(data) %like% ".y"])]
data <- data[complete.cases(data), ]
colnames(data) <- gsub(".y", "", colnames(data))
data$sample <- gsub("(_WG).*", "\\1", data$sample)

Data <- aggregate(data[, 3:ncol(data)], list(data$sample, data$ActualClass), mean)
Data$sample <- Data$Group.1
y <- Data$Group.2
y <- factor(y, levels = c("negative", "positive"))
Data <- Data[, !(colnames(Data) %in% c("Group.1", "Group.2"))]

Data_butler <- aggregate(butler[, 3:ncol(butler)], list(butler$sample, butler$ActualClass), mean)
row.names(Data_butler) <- Data_butler$Group.1
Data_butler <- Data_butler[, names]

### Run classifier
algorithm <- "rf"
source(class)
outcomes_lfs <- All.kFold
outcomes_neg_butler <- But.kfold
saveRDS(outcomes_lfs, file.path(path, "outcomes_integrated_status.rds"))
saveRDS(outcomes_neg_butler, file.path(path, "outcomes_integrated_status_butler.rds"))
