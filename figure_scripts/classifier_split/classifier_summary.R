library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(ggpubr)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split"

types <- c("status_test", "status_val", "lfs_test", "lfs_val")
titles <- c("Validation Cohort", "Test Cohort", "Validation Cohort", "Test Cohort")
names <- c("Healthy Non-Carriers vs LFS Carriers\nValidation Cohort", "Healthy Non-Carriers vs LFS Carriers\nTest Cohort",
           "LFS Cancer Free vs Active Cancer\nValidation Cohort", "LFS Cancer Free vs Active Cancer\nTest Cohort")

samples <- read.delim("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split/samples_split.txt")
samples_status <- samples[samples$cancer_status %in% c("negative", "healthy"), ]
samples_status$ActualClass <- ifelse(samples_status$cancer_status == "negative", "positive", "negative")

samples_lfs <- samples[samples$cancer_status %in% c("negative", "positive"), ]
samples_lfs$ActualClass <- samples_lfs$cancer_status

data_comp <- data.frame()

for (i in c(1:length(types))) {
  ### Set variables
  type <- types[[i]]
  name <- names[[i]]
  title <- titles[[i]]
  
  ### Find files
  filenames <- list.files(path, paste0(type, ".rds"), full.names = TRUE)
  #filenames <- filenames[!(grepl("integrated", filenames))]
  analyses <- list.files(path, paste0(type, ".rds"), full.names = FALSE)
  #analyses <- analyses[!(grepl("integrated", analyses))]
  analyses <- gsub("outcomes_", "", analyses)
  analyses <- sub(paste0("_", type, ".rds"), "", analyses)
  
  data_roc <- data.frame()
  data_auc <- data.frame()
  data_cm <- data.frame()
  
  for (j in c(1:length(filenames))) {
    ### Set variables
    file <- filenames[[j]]
    analysis <- analyses[[j]]
    
    ### Read in file
    data <- readRDS(file)
    
    ### Calculuate confidence intervals
    if(type %like% "test") {
    ci <- c()
    for (k in 1:length(data)) {
      x <- unlist(data[k], recursive = FALSE)
      x <- bind_rows(x)
      y <- auc(x$ActualClass, x$positive)
      ci <- c(ci, y)
    }
    auc <- format(signif(mean(ci), digits = 3), nsmall = 3)
    se <- sd(ci)/sqrt(length(ci))
    qt <- qt((0.95+1)/2, df = length(ci) - 1)
    lower <- format(signif(mean(ci) - qt*se, digits = 3), nsmall = 3)
    upper <- format(signif(mean(ci) + qt*se, digits = 3), nsmall = 3)
    }
    
    if(type == "status_val") {
      ci <- c()
      for (k in 1:length(data)) {
        x <- unlist(data[k], recursive = FALSE)
        x <- bind_rows(x)
        x <- merge(x, samples_status[, c("sWGS", "ActualClass")], by.x = "sample", by.y = "sWGS")
        y <- auc(x$ActualClass, x$positive)
        ci <- c(ci, y)
      }
      auc <- format(signif(mean(ci), digits = 3), nsmall = 3)
      se <- sd(ci)/sqrt(length(ci))
      qt <- qt((0.95+1)/2, df = length(ci) - 1)
      lower <- format(signif(mean(ci) - qt*se, digits = 3), nsmall = 3)
      upper <- format(signif(mean(ci) + qt*se, digits = 3), nsmall = 3)
    }
    
    if(type == "lfs_val") {
      ci <- c()
      for (k in 1:length(data)) {
        x <- unlist(data[k], recursive = FALSE)
        x <- bind_rows(x)
        x <- merge(x, samples_lfs[, c("sWGS", "ActualClass")], by.x = "sample", by.y = "sWGS")
        y <- auc(x$ActualClass, x$positive)
        ci <- c(ci, y)
      }
      auc <- format(signif(mean(ci), digits = 3), nsmall = 3)
      se <- sd(ci)/sqrt(length(ci))
      qt <- qt((0.95+1)/2, df = length(ci) - 1)
      lower <- format(signif(mean(ci) - qt*se, digits = 3), nsmall = 3)
      upper <- format(signif(mean(ci) + qt*se, digits = 3), nsmall = 3)
    }
    
    ### Unlist files
    data <- unlist(data, recursive = FALSE)
    data <- bind_rows(data)
    
    if(type == "status_val") {
      data <- merge(data, samples_status[, c("sWGS", "ActualClass")], by.x = "sample", by.y = "sWGS")
      #data <- aggregate(data[, c("positive")], list(data$sample, data$ActualClass), mean)
      #colnames(data) <- c("sample", "ActualClass", "positive")
      #data$negative <- 1-data$positive
      #data$PredictedClass <- ifelse(data$positive >=0.5, "positive", "negative")
    }
    
    if(type == "lfs_val") {
      data <- merge(data, samples_lfs[, c("sWGS", "ActualClass")], by.x = "sample", by.y = "sWGS")
      #data <- aggregate(data[, c("positive")], list(data$sample, data$ActualClass), mean)
      #colnames(data) <- c("sample", "ActualClass", "positive")
      #data$negative <- 1-data$positive
      #data$PredictedClass <- ifelse(data$positive >=0.5, "positive", "negative")
    }

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
    
    #auc <- auc(data$ActualClass, data$positive)
    #auc <- format(signif(auc, digits = 3), nsmall = 3)
    auc <- data.frame(Metric = analysis, 
                      auc = auc,
                      upper = upper,
                      lower = lower)
    
    ### Append to dataframes
    data_roc <- rbind(data_roc, roc)
    data_auc <- rbind(data_auc, auc)
    data_cm <- rbind(data_cm, cm)
  }
  ### Format plotting tables
  labels <- c("Breakpoint", "Dinucleotide", "End Motif", "Fragment Length", "Integrated", "Nucleosome Footprint", 
              "Fragment Ratio", "Accessibility TCGA/DHS", "Accessibility TFBS")
  data_roc$Metric <- factor(data_roc$Metric, levels = unique(data_roc$Metric),
                            labels = labels)
  data_auc$Metric <- factor(data_auc$Metric, levels = unique(data_auc$Metric),
                            labels = labels)
  data_cm$Metric <- factor(data_cm$Metric, levels = unique(data_cm$Metric),
                           labels = labels)
  data_auc <- data_auc[order(data_auc$auc), ]
  order <- data_auc$Metric
  data_auc$Metric <- factor(data_auc$Metric, levels = order)
  data_auc$comparison <- type
  data_auc$cohort <- title
  
  data_cm$Metric <- factor(data_cm$Metric, levels = order)
  
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
    geom_text(data = data_auc, aes(x = 1, y = c(seq(0, 0.49, 0.06)), color = Metric, 
                                   label = paste0(Metric, " = ", auc, " (", lower, ",", upper, ")")), 
              size = 4, hjust = 1) +
    scale_color_manual(values = c("Integrated" = "#000000", "End Motif" = "#009E73", "Nucleosome Footprint" = "#E69F00", "Fragment Ratio" = "#0072B2", "Breakpoint" = "grey65",
                                  "Fragment Length" = "#D55E00", "Accessibility TFBS" = "#CC79A7", "Dinucleotide" = "#6A3D9A", "Accessibility TCGA/DHS" = "#56B4E9")) +
    xlab("False Positive Rate") + 
    ylab("True Positive Rate") +
    ggtitle(title) + 
    theme +
    theme(legend.position = "none")
  AUC_plot
  
  ### Plot Performance Metrics
  CM_plot <- ggplot(data_cm) +
    geom_point(aes(sensitivity, specificity, color = Metric), size = 3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
    scale_color_manual(values = c("Integrated" = "#000000", "End Motif" = "#009E73", "Nucleosome Footprint" = "#E69F00", "Fragment Ratio" = "#0072B2", "Breakpoint" = "grey65",
                                  "Fragment Length" = "#D55E00", "Accessibility TFBS" = "#CC79A7", "Dinucleotide" = "#6A3D9A", "Accessibility TCGA/DHS" = "#56B4E9")) +
    xlab("Sensitivity") + 
    ylab("Specificity") +
    labs(color = "") +
    ggtitle(title) + 
    theme +
    theme(legend.position = c(0.5, 0.21)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    guides(color = guide_legend(ncol = 2))
  CM_plot
  
  ### Save Figure
  #fig <- ggarrange(AUC_plot, CM_plot, align = "h")
  #annotate_figure(fig, top = text_grob(name, face = "bold", size = 14))
  ggsave(file.path(path, paste0("classifier_auc_", type, ".pdf")), AUC_plot,  width = 5, height = 4.5)
  ggsave(file.path(path, paste0("classifier_cm_", type, ".pdf")), CM_plot,  width = 4.5, height = 4.5)
  
  ### Apphend comparison data
  data_comp <- rbind(data_comp, data_auc)
}

### Compare AUCs between test and validation sets
data_comp$cohort <- gsub(" Cohort", "", data_comp$cohort)
data_comp$auc <- as.numeric(data_comp$auc)

data_comp1 <- data_comp[c(1:(nrow(data_comp)/2)), ]
data_comp1 <- data_comp1[order(data_comp1$Metric), ]
data_comp1_stats <- data_comp1 %>%
  group_by(Metric) %>%
  mutate(diff = auc - lag(auc, default = first(auc)))
data_comp1$cohort <- factor(data_comp1$cohort, levels = c("Validation", "Test"))

data_comp2 <- data_comp[(nrow(data_comp)/2 + 1):nrow(data_comp), ]
data_comp2 <- data_comp2[order(data_comp2$Metric), ]
data_comp2_stats <- data_comp2 %>%
  group_by(Metric) %>%
  mutate(diff = auc - lag(auc, default = first(auc)))
data_comp2$cohort <- factor(data_comp2$cohort, levels = c("Validation", "Test"))

comp_plot1 <- ggplot(data_comp1) +
  geom_point(aes(cohort, auc, color = Metric), size = 3) +
  geom_line(aes(cohort, auc, group = Metric, color = Metric), size = 1) +
  scale_color_manual(values = c("Integrated" = "#000000", "End Motif" = "#009E73", "Nucleosome Footprint" = "#E69F00", "Fragment Ratio" = "#0072B2", "Breakpoint" = "grey65",
                                "Fragment Length" = "#D55E00", "Accessibility TFBS" = "#CC79A7", "Dinucleotide" = "#6A3D9A", "Accessibility TCGA/DHS" = "#56B4E9")) +
  xlab("Cohort") + 
  ylab("AUC-ROC") +
  labs(color = "") +
  ggtitle("AUC Comparison") + 
  theme +
  theme(legend.position = c(0.5, 0.21)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(expand = expansion(mult=c(0.25,0.25))) +
  guides(color = guide_legend(ncol = 2))
comp_plot1

ggsave(file.path(path, paste0("classifier_auc_comparison_status.pdf")), comp_plot1,  width = 4.5, height = 4.5)


comp_plot2 <- ggplot(data_comp2) +
  geom_point(aes(cohort, auc, color = Metric), size = 3) +
  geom_line(aes(cohort, auc, group = Metric, color = Metric), size = 1) +
  scale_color_manual(values = c("Integrated" = "#000000", "End Motif" = "#009E73", "Nucleosome Footprint" = "#E69F00", "Fragment Ratio" = "#0072B2", "Breakpoint" = "grey65",
                                "Fragment Length" = "#D55E00", "Accessibility TFBS" = "#CC79A7", "Dinucleotide" = "#6A3D9A", "Accessibility TCGA/DHS" = "#56B4E9")) +
  xlab("Cohort") + 
  ylab("AUC-ROC") +
  labs(color = "") +
  ggtitle("AUC Comparison") + 
  theme +
  theme(legend.position = c(0.5, 0.21)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(expand = expansion(mult=c(0.25,0.25))) +
  guides(color = guide_legend(ncol = 2))
comp_plot2

ggsave(file.path(path, paste0("classifier_auc_comparison_lfs.pdf")), comp_plot2,  width = 4.5, height = 4.5)
