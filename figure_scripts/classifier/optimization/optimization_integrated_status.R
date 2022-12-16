library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(multiROC)

### Set variables
path <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization"

### List classifiers
class_list <- list.files(path, "status.rds", full.names = TRUE)
class_list <- class_list[!grepl("integrated", class_list)]
class_list <- class_list[!grepl("end_motif", class_list)]

names <- list.files(path, "status.rds", full.names = FALSE)
names <- names[!grepl("integrated", names)]
names <- names[!grepl("end_motif", names)]
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

data <- data[complete.cases(data), ]
y <- data$ActualClass
data <- data[, !(colnames(data) %in% c("ActualClass", "sample"))]

### Set models to run
#Features.CVparam <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = TRUE, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE)
classifiers <- c("rf", "glm", "gbm", "svmRadial", "knn")
names <- c("RF", "GLM", "GBM", "SVM", "KNN")

### Run tests (healthy vs LFS)
results <- data.frame(test = c(), kappa = c(), CI = c())

for (k in c(1:length(classifiers))){
  ### Make table
  Data <- data
  
  y <- factor(y, levels = c("negative", "positive"))
  
  ### Set variables
  class <- classifiers[[k]]
  name <- names[[k]]
  
  ### Run suite of classifiers
  set.seed(123)
  source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/classifier_opt.R")
  
  ### Get Kappa performance
  mean <- mean(kappa)
  margin <- qt(0.975, df = length(kappa) - 1)*sd(kappa)/sqrt(length(kappa))
  
  ### Apphend to dataframe
  a <- data.frame(test = name, kappa = mean, error = margin)
  results <- rbind(results, a)
}

### Set theme
scaleFUN <- function(x) sprintf("%.2f", x)
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.position = "none",
               legend.background = element_blank(),
               axis.text = element_text(size = 12),
               axis.title = element_text(size = 12),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Plot Kappa values
results <- results[order(results$kappa), ]
results$test <- factor(results$test, levels = results$test)
plot <- ggplot(results) +
  geom_point(aes(x = test, y = kappa)) +
  geom_errorbar(aes(x = test, ymin = kappa - error, ymax = kappa + error), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  xlab("Algorithm") + 
  ylab("Kappa") +
  labs(color = "") + 
  ggtitle("Integrated Healthy Control\nvs LFS Cancer Free") + 
  theme +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip()
plot
ggsave(file.path(outdir, paste0("kappa_integrated_status.pdf")), plot, width = 3, height = 3.5)
