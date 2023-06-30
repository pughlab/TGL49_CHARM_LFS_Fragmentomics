library(caret)
library(dplyr)
library(tidyverse)
library(pROC)

### Set variables
source("figure_scripts/classifier/FuncClassifier.R")
class <- "figure_scripts/classifier/classifier.R"
path <- "data/end_motifs"
healthy_path <- "hbc/end_motifs"
outdir <- ""

### Import data 
data_end <- read.delim(list.files(path, "CHARM_LFS_genome2_end_motifs.txt", full.names = TRUE))
normal_end <- read.delim(list.files(healthy_path, "CHARM_HBC_genome2_end_motifs.txt", full.names = TRUE))
data_samples <- read.delim("sample_list.txt")

group_lfs <- read.delim(list.files("data/read_group", ".txt", full.names = TRUE))
group_hbc <- read.delim(list.files("hbc/read_group", ".txt", full.names = TRUE))

### Remove failed data_samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_end), ]

### Find runs with both LFS and HBCs
group_lfs$type <- "LFS"
group_hbc$type <- "HBC"
groups <- bind_rows(group_lfs, group_hbc)
groups <- groups %>%
  group_by(type, V1) %>%
  dplyr::summarise()
groups <- groups[groups$V1 %in% groups$V1[duplicated(groups$V1)], ]

group_lfs <- group_lfs[group_lfs$V1 %in% groups$V1, ]
group_hbc <- group_hbc[group_hbc$V1 %in% groups$V1, ]

### Set sample lists
samples_neg <- data_samples[data_samples$cancer_status == "negative", ]
samples_neg <- data_samples[data_samples$sWGS %in% group_lfs$sample, ]

### Format table for classifier
data_end[is.na(data_end)] <- 0
row.names(data_end) <- data_end$motif
data_end <- data_end[, colnames(data_end) %in% data_samples$sWGS]
data_end <- as.data.frame(t(data_end))
data_end <- data_end[samples_neg$sWGS, ]

normal_end[is.na(normal_end)] <- 0
row.names(normal_end) <- normal_end$motif
normal_end <- normal_end[, -1]
normal_end <- as.data.frame(t(normal_end))
normal_end <- normal_end[row.names(normal_end) %in% group_hbc$sample, ]

### Make Healthy vs LFS Negative
Data <- data_end
Data <- bind_rows(Data, normal_end)
Data$sample <- row.names(Data)

y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
y <- factor(y, levels = c("negative", "positive"))

algorithm <- "gbm"

### Run the classifier
Features.CVparam <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = TRUE, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE)

gbmGrid <- expand.grid(interaction.depth = 3,
                       n.trees = 150,
                       shrinkage = 0.1,
                       n.minobsinnode = 2)

# Split the data into 10 folds
Splits10 <- list()
All.kFold <- list()

for (j in 1:10) {
  # split the data for 100 times, result of each split is saved in Splits10
  Splits10[[j]] <- SplitkFold(Data, y, 10)
  
  #10-fold cross validation to calculate the probability (score) of each sample being cancer
  kFold.list <- list()
  kFold.but <- list()
  
  for(i in 1:10) {
    Indices = Splits10[[j]]$samples[[i]]
    classes.df = Splits10[[j]]$df
    
    TrainData <- Splits10[[j]][["data"]][Indices, ]
    TrainPheno <- classes.df[Indices,]
    
    TestData <- Splits10[[j]][["data"]][!(row.names(Splits10[[j]][["data"]]) %in% row.names(TrainData)), ]
    TestPheno <- classes.df[classes.df$ID %in% row.names(TestData), ]
    
    Model <- train(x = TrainData, y = TrainPheno$Classes, trControl = Features.CVparam, method = algorithm, metric = "Kappa", tuneGrid = gbmGrid)
    Prediction.classProbs <- predict(Model, newdata = TestData, type = "prob") %>% data.frame
    
    Prediction.classProbs$ActualClass <- TestPheno$Classes
    Prediction.classProbs$PredictedClass <- predict(Model, newdata = TestData, type = "raw")
    Prediction.classProbs$sample <- row.names(TestData)
    
    kFold.list[[i]] <- Prediction.classProbs
  }
  
  Predicted.classProbs <- kFold.list[[1]]$TestPred
  for (i in 2:10) {
    Predicted.classProbs <- bind_rows(Predicted.classProbs, kFold.list[[i]]$TestPred)
  }
  
  All.kFold[[j]] <- kFold.list
}

outcomes <- All.kFold

### Find AUC and Performance
### Unlist files
outcomes <- unlist(outcomes, recursive = FALSE)
outcomes <- bind_rows(outcomes)

### Make confusion matrix
outcomes$ActualClass <- factor(outcomes$ActualClass, levels = c("negative", "positive"))
outcomes$PredictedClass <- factor(outcomes$PredictedClass, levels = c("negative", "positive"))
cm <- confusionMatrix(outcomes$PredictedClass, outcomes$ActualClass)
cm <- data.frame(Metric = "End Motif",
                 sensitivity = cm[["byClass"]][["Sensitivity"]],
                 specificity = cm[["byClass"]][["Specificity"]])
outcomes$ActualClass <- ifelse(outcomes$ActualClass == "negative", 0, 1)

roc <- roc(outcomes$ActualClass, outcomes$positive)
specificity <- roc[["specificities"]]
sensitivity <- roc[["sensitivities"]]
sensitivity <- 1 - sensitivity

roc <- as.data.frame(cbind(sensitivity, specificity))
colnames(roc) <- c("sensitivity", "specificity")
roc$Metric <- analysis

auc <- auc(outcomes$ActualClass, outcomes$positive)
auc <- format(signif(auc, digits = 3), nsmall = 3)

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.title = element_text(size = 13),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.position = "right",
               legend.background = element_blank(),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

### Plot the ROCs
AUC_plot <- ggplot() +
  geom_line(data = roc, aes(sensitivity, specificity)) +
  geom_abline(intercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_text(aes(x = 0.7, y = 0.1, label = paste0("AUC = ", auc)), size = 4) +
  scale_color_manual(values = c("grey65", "#000000", "#009E73", "#E69F00", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9")) +
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") +
  ggtitle("LFS vs HBC samples on same run") + 
  theme +
  theme(legend.position = "none")
AUC_plot

ggsave(file.path(outdir, "classifier_end_motif_investigation.pdf"), AUC_plot, width = 4.5, height = 4.5)

