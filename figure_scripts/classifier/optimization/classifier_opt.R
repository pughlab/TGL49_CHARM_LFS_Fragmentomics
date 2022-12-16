library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(NMF)
library(Rtsne)
library(caret)
library(randomForest)
library(glmnet)
library(limma)

############################# load the input data here ###############################
# Data = n_sample * n_feature dataframe
# y = vector of length n_sample, each element represents the label of the sample (positive or negative)

############################# classifer starts here###############################

Features.CVparam <- trainControl(method = "repeatedcv", number = 10, repeats = 1, verboseIter = TRUE, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE)

# Split the data into 10 folds
Splits <- list()
kappa <- c()

for (i in 1:10) {
  # split the data for 100 times, result of each split is saved in Splits10
  Splits[[i]] <- downSample(Data, y, list = TRUE)
}

for(j in 1:length(Splits)) {

    TrainData <- Splits[[j]][["x"]]
    TrainPheno <- Splits[[j]][["y"]]

    Model <- train(x = TrainData, y = TrainPheno, trControl = Features.CVparam, method = class, metric = "Kappa")

    kappa <- c(kappa, Model[["results"]][["Kappa"]])
}


# Features.list <- list()
#
# for (j in 1:100) {
#   # Get features from the 100 repeats, each run has 10 folds (runs)
#   Feature.onetime <- c()
#   kFold.list <- All.kFold[[j]]
#   for (i in 1:10) {
#     Feature.onetime <- c(Feature.onetime, kFold.list[[i]]$Model$finalModel$xNames)
#   }
#   Feature.onetime <- unique(Feature.onetime)
#
#   Features.list[[j]] <- Feature.onetime
# }
#
# saveRDS(Features.list, 'Features.list.rds') # 100 groups of features saved here
#
# Features <- unique(unlist(Features.list))
# write.table(Features, 'Features.tsv', quote=F, col.names=F, row.names=F)
