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

Features.CVparam <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = TRUE, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE)

# Split the data into 10 folds
Splits10 <- list()
All.kFold <- list()
But.kfold <- list()

for (j in 1:100) {
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

    Model <- train(x = TrainData, y = TrainPheno$Classes, trControl = Features.CVparam, method = algorithm, metric = "Kappa")
    Prediction.classProbs <- predict(Model, newdata = TestData, type = "prob") %>% data.frame

    Prediction.classProbs$ActualClass <- TestPheno$Classes
    Prediction.classProbs$PredictedClass <- predict(Model, newdata = TestData, type = "raw")
    Prediction.classProbs$sample <- row.names(TestData)

    kFold.list[[i]] <- Prediction.classProbs
    
    #Predict Butler data
    Prediction.but <- predict(Model, newdata = Data_butler, type = "prob") %>% data.frame
    
    Prediction.but$PredictedClass <- predict(Model, newdata = Data_butler, type = "raw")
    Prediction.but$sample <- row.names(Data_butler)
    
    kFold.but[[i]] <- Prediction.but
  }

  Predicted.classProbs <- kFold.list[[1]]$TestPred
  for (i in 2:10) {
    Predicted.classProbs <- bind_rows(Predicted.classProbs, kFold.list[[i]]$TestPred)
  }

  All.kFold[[j]] <- kFold.list
  But.kfold[[j]] <- kFold.but
}
