library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(ggpubr)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

types <- c("status", "lfs")
names <- c("Healthy Non-Carriers vs LFS Carriers", "LFS Cancer Free vs Active Cancer")

for (i in c(1:length(types))) {
  ### Set variables
  type <- types[[i]]
  name <- names[[i]]
  
  ### Find files
  filenames <- list.files(path, paste0(type, ".rds"), full.names = TRUE)
  filenames <- filenames[!(grepl("integrated", filenames))]
  analyses <- list.files(path, paste0(type, ".rds"), full.names = FALSE)
  analyses <- analyses[!(grepl("integrated", analyses))]
  analyses <- gsub("outcomes_", "", analyses)
  analyses <- sub("_[^_]+$", "", analyses)
  
  data_roc <- data.frame()
  data_auc <- data.frame()
  data_cm <- data.frame()
  
  for (j in c(1:length(filenames))) {
    ### Set variables
    file <- filenames[[j]]
    analysis <- analyses[[j]]
    
    ### Read in file
    data <- readRDS(file)
    
    ### Unlist files
    data <- unlist(data, recursive = FALSE)
    data <- bind_rows(data)

    ### Make confusion matrix
    data$ActualClass <- factor(data$ActualClass, levels = c("negative", "positive"))
    data$PredictedClass <- factor(data$PredictedClass, levels = c("negative", "positive"))
    cm <- confusionMatrix(data$PredictedClass, data$ActualClass)
    cm <- data.frame(Metric = analysis,
                     sensitivity = cm[["byClass"]][["Sensitivity"]],
                     specificity = cm[["byClass"]][["Specificity"]])
    
    ### Format table
    data$ActualClass <- ifelse(data$ActualClass == "negative", 0, 1)
    
    ### Calculate auc and roc values
    roc <- roc(data$ActualClass, data$positive)
    specificity <- roc[["specificities"]]
    sensitivity <- roc[["sensitivities"]]
    sensitivity <- 1 - sensitivity
    
    roc <- as.data.frame(cbind(sensitivity, specificity))
    colnames(roc) <- c("sensitivity", "specificity")
    roc$Metric <- analysis
    
    auc <- auc(data$ActualClass, data$positive)
    auc <- format(signif(auc, digits = 3), nsmall = 3)
    auc <- data.frame(Metric = analysis, auc = auc)
    
    ### Append to dataframes
    data_roc <- rbind(data_roc, roc)
    data_auc <- rbind(data_auc, auc)
    data_cm <- rbind(data_cm, cm)
  }
  ### Format plotting tables
  data_roc$Metric <- factor(data_roc$Metric, levels = unique(data_roc$Metric),
                            labels = c("Breakpoint", "Dinucleotide", "End Motif", "Fragment Length", "Nucleosome Distance", 
                                       "Fragment Ratio", "Accessibility TCGA/DHS", "Accessibility TFBS"))
  data_auc$Metric <- factor(data_auc$Metric, levels = unique(data_auc$Metric),
                            labels = c("Breakpoint", "Dinucleotide", "End Motif", "Fragment Length", "Nucleosome Distance", 
                                       "Fragment Ratio", "Accessibility TCGA/DHS", "Accessibility TFBS"))
  data_cm$Metric <- factor(data_cm$Metric, levels = unique(data_cm$Metric),
                           labels = c("Breakpoint", "Dinucleotide", "End Motif", "Fragment Length", "Nucleosome Distance", 
                                      "Fragment Ratio", "Accessibility TCGA/DHS", "Accessibility TFBS"))
  
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
  AUC_plot <- ggplot(data_roc) +
    geom_line(aes(sensitivity, specificity, group = Metric, color = Metric)) +
    geom_abline(intercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_text(data = data_auc, aes(x = 0.7, y = c(seq(0, 0.42, 0.06)), color = Metric, label = paste0(Metric, " = ", auc)), size = 4) +
    scale_color_manual(values = c("grey65", "#000000", "#009E73", "#E69F00", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9")) +
    xlab("False Positive Rate") + 
    ylab("True Positive Rate") +
    ggtitle("ROC Curves") + 
    theme +
    theme(legend.position = "none")
  AUC_plot
  
  ### Plot Performance Metrics
  CM_plot <- ggplot(data_cm) +
    geom_point(aes(sensitivity, specificity, color = Metric), size = 3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
    scale_color_manual(values = c("grey65", "#000000", "#009E73", "#E69F00", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9")) +
    xlab("Sensitivity") + 
    ylab("Specificity") +
    labs(color = "") +
    ggtitle("Performance") + 
    theme +
    theme(legend.position = c(0.5, 0.21)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    guides(color = guide_legend(ncol = 2))
  CM_plot
  
  ### Save Figure
  fig <- ggarrange(AUC_plot, CM_plot, align = "h")
  annotate_figure(fig, top = text_grob(name, face = "bold", size = 14))
  ggsave(file.path(path, paste0("classifier_performance_", type, ".pdf")), width = 9, height = 4.5)
}

### Make Integrated ROCs




