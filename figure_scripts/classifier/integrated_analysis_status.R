library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(ggpubr)
library(ComplexHeatmap)
library(rstatix)

### Set variables
path <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Find files
filenames <- list.files(path, "status.rds", full.names = TRUE)
filenames <- filenames[grepl("integrated", filenames)]

### Read in file
data <- readRDS(filenames)
samples <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Unlist files
data <- unlist(data, recursive = FALSE)
data <- bind_rows(data)

### Make confusion matrix
data$ActualClass <- factor(data$ActualClass, levels = c("negative", "positive"))
data$PredictedClass <- factor(data$PredictedClass, levels = c("negative", "positive"))
cm <- confusionMatrix(data$PredictedClass, data$ActualClass)
cm <- data.frame(sensitivity = cm[["byClass"]][["Sensitivity"]],
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

auc <- auc(data$ActualClass, data$positive)
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
               legend.text = element_text(size = 11),
               legend.position = "none",
               legend.background = element_blank(),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

### Plot the ROCs
AUC_plot <- ggplot(roc) +
  geom_line(aes(sensitivity, specificity)) +
  geom_abline(intercept = 0, linetype = "dashed", alpha = 0.5) +
  annotate(geom = "text", x = 0.7, y = 0.06, label = paste0("Integrated AUC = ", auc), size = 4) +
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") +
  ggtitle("Integrated Healthy Control\nvs LFS Cancer Free") + 
  theme
AUC_plot

ggsave(file.path(path, "classifier_performance_integrated_status.pdf"), AUC_plot, width = 4.5, height = 4.5)

### Read in other metrics and combine
filenames <- list.files(path,  "status.rds", full.names = TRUE)
filenames <- filenames[!(grepl("integrated", filenames))]
filenames <- filenames[!(grepl("ratio", filenames))]
filenames <- filenames[!(grepl("end_motif", filenames))]
analyses <- list.files(path, "status.rds", full.names = FALSE)
analyses <- analyses[!(grepl("integrated", analyses))]
analyses <- analyses[!(grepl("ratio", analyses))]
analyses <- analyses[!(grepl("end_motif", analyses))]
analyses <- gsub("outcomes_", "", analyses)
analyses <- sub("_[^_]+$", "", analyses)

data_scores <- data[, c("positive", "ActualClass", "sample")]
data_scores$ActualClass <- factor(data_scores$ActualClass, levels = c(0,1),
                                  labels = c("negative", "positive"))
data_scores <- data_scores %>%
  group_by(sample) %>%
  dplyr::summarise(score=mean(positive),
                   ActualClass=unique(ActualClass))

for (i in c(1:length(filenames))) {
  ### Set variables
  file <- filenames[[i]]
  analysis <- analyses[[i]]
  
  ### Read in file
  data <- readRDS(file)
  
  ### Unlist files
  data <- unlist(data, recursive = FALSE)
  data <- bind_rows(data)
  
  ### Summarize outputs
  data$ActualClass <- factor(data$ActualClass, levels = c("negative", "positive"))
  data <- data %>%
    group_by(sample) %>%
    dplyr::summarise(score=mean(positive),
                     ActualClass=unique(ActualClass))
  
  ### Merge to data
  data_scores <- merge(data_scores, data, by = c("ActualClass", "sample"), all = TRUE)
}

data_scores <- data_scores[complete.cases(data_scores), ]
colnames(data_scores) <- c("ActualClass", "sample", "integrated", analyses)

### Merge with clinical data and format
samples <- samples[samples$sWGS %in% data_scores$sample, ]
data <- merge(data_scores, samples[, c("ext_ID", "timepoint", "family", "Age", "previous_cancer", "germline_mutation", "mutation_type", "sWGS")], by.x = "sample", by.y = "sWGS", all = TRUE)
data[is.na(data)] <- "Healthy"

data$ActualClass <- factor(data$ActualClass, levels = c("negative", "positive"),
                           labels = c("Healthy", "LFS"))
data$Age <- factor(data$Age,levels = c("Healthy", "adult", "pediatric"),
                       labels = c("Healthy", "Adult", "Pediatric"))
data$previous_cancer <- factor(data$previous_cancer, levels = c("Healthy", "yes", "no", ""),
                                   labels = c("Healthy", "Yes", "No", "Unknown"))
data$mutation_type <- factor(data$mutation_type, levels = c("Healthy", "LOF", "Splice", "1", "2", "3", "5"),
                                 labels = c("Healthy", "LOF", "Splice", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 5"))
data <- data[complete.cases(data), ]

### Make matrixes
data_class <- as.matrix(t(data[,c(4:8)]))

### Set colours
col_diag <- c(LFS = "grey95", Healthy = "#B2DF8A")
col_previous <- c(Healthy = "#B2DF8A", Yes = "#fb9a99", No = "grey95")
col_age <- c(Healthy = "#B2DF8A", Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_type <- c(Healthy = "#B2DF8A", LOF = "#FC8D62", Splice = "#66C2A5", "Cluster 1" = "#E5C494", "Cluster 2" = "#FFD92F", "Cluster 3" = "#8DA0CB", "Cluster 5" = "#E78AC3")

### Set annotations
top_annotation <- HeatmapAnnotation("Metrics" = anno_barplot(t(data_class),
                                                             axis_param = list(gp = gpar(fontsize = 6)),
                                                             border = FALSE),
                                    "Integrated" = anno_points(data$integrated),
                                    border = FALSE,
                                    height = unit(1.75, "cm"),
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 8),
                                    simple_anno_size = unit(0.25, "cm"))

bot_annotation <- HeatmapAnnotation("Germline Mutation" = data$mutation_type,
                                    "Cancer History" = data$previous_cancer,
                                    "Patient Type" = data$Age,
                                    "Diagnosis" = data$ActualClass,
                                    col = list("Germline Mutation" = col_type,
                                               "Cancer History" = col_previous,
                                               "Patient Type" = col_age,
                                               "Diagnosis" = col_diag),
                                    border = FALSE,
                                    height = unit(1.75, "cm"),
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "right",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 8),
                                    simple_anno_size = unit(0.25, "cm"))

### Make heatmap
Heatmap <- Heatmap(data_class,
                   #col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   top_annotation = top_annotation,
                   bottom_annotation = bot_annotation,
                   row_labels = NULL,
                   column_title_gp = gpar(fontsize = 6),
                   row_title_gp = gpar(fontsize = 0),
                   row_title_rot = 0,
                   row_dend_width = unit(0, "cm"),
                   border = FALSE,
                   height = unit(0, "cm"))
Heatmap

### Calculate statistics
data_stats <- data %>%
  group_by(mutation_type) %>%
  dplyr::summarise(mean=mean(integrated),
                   sd=sd(integrated),
                   N=n(),
                   patients=length(unique(ext_ID)))
data_stats$annot <- c("", rep("***", nrow(data_stats) - 1))

### Plot
data <- data %>%
  group_by(mutation_type) %>%
  arrange(integrated)
data_stats$Xn <- c(1:nrow(data_stats))
fig_type <- ggplot(data) +
  stat_summary(aes(mutation_type, integrated), 
               fun = median, geom = "crossbar", size = 0.2, width = 0.5, group = 1, show.legend = FALSE, color = "red") +
  geom_rect(data = data_stats, aes(xmin = Xn - 0.5, xmax = Xn + 0.5, ymin = -Inf, ymax = Inf, fill = mutation_type), alpha = 0.5, stat = "identity") +
  geom_point(aes(mutation_type, integrated), alpha = 0.75,
             position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats, aes(x = mutation_type, y = -0.15, label = N), size = 4) +
  geom_text(data = data_stats, aes(x = mutation_type, y = -0.05, label = paste0("(", patients, ")")), size = 4) +
  geom_text(data = data_stats, aes(x = mutation_type, y = 1.1, label = annot), size = 5) +
  ggtitle("Germline Mutation Type") + 
  xlab("") +
  ylab("Integrated Score") +
  scale_fill_manual(values = c("white", "#FB9A99", "#A6CEE3", "#FDBF6F","#CAB2D6", "grey65", "#B2DF8A")) +
  scale_color_manual(values = c("black", "red")) +
  labs(fill = "") +
  theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(limits=c(-0.2, 1.15), expand = c(0,0))
fig_type

