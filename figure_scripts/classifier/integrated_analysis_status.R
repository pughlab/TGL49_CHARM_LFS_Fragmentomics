library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(ggpubr)
library(ComplexHeatmap)
library(rstatix)
library(circlize)

### Set variables
path <- "data/classifier"

### Find files
filenames <- list.files(path, "status.rds", full.names = TRUE)
filenames <- filenames[grepl("integrated", filenames)]

### Read in file
data <- readRDS(filenames)
samples <- read.delim("sample_list.txt")
samples_healthy <- read.delim("hbc/sample_list.txt")

cancer_patients <- unique(samples$ext_ID[samples$cancer_status == "positive"])

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

### Read in other metrics and combine
filenames <- list.files(path,  "status.rds", full.names = TRUE)
filenames <- filenames[!(grepl("integrated", filenames))]
#filenames <- filenames[!(grepl("ratio", filenames))]
filenames <- filenames[!(grepl("end_motif", filenames))]
analyses <- list.files(path, "status.rds", full.names = FALSE)
analyses <- analyses[!(grepl("integrated", analyses))]
#analyses <- analyses[!(grepl("ratio", analyses))]
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
data <- merge(data_scores, samples[, c("ext_ID", "timepoint", "family", "Age", "years", "previous_cancer", "germline_mutation", "mutation_type", "sWGS")], 
              by.x = "sample", by.y = "sWGS", all = TRUE)
data[is.na(data)] <- "Healthy"

data$ActualClass <- factor(data$ActualClass, levels = c("negative", "positive"),
                           labels = c("Healthy", "LFS"))
data$Age <- factor(data$Age,levels = c("Healthy", "adult", "pediatric"),
                       labels = c("Healthy", "Adult", "Pediatric"))
data$previous_cancer <- factor(data$previous_cancer, levels = c("Healthy", "yes", "no", ""),
                                   labels = c("Healthy", "Yes", "No", "Unknown"))
data$mutation_type <- factor(data$mutation_type, levels = c("Healthy", "LOF", "Splice", "1", "2", "3", "5"),
                                 labels = c("Healthy", "LOF", "Splice Site", "Missense 1", "Missense 2", "Missense 3", "Missense 5"))
data <- data[complete.cases(data), ]
data <- data[order(data$integrated), ]

### Make heatmap matrix
data_class <- as.matrix(t(data[,c(4:10)]))
row.names(data_class) <- c("Breakpoint", "Dinucleotide", "Fragment Length", "Nucleosome Footprint", "Fragment Ratio", "Accessibility TCGA/DHS", "Accessibility TFBS")

### Set colours
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
col_diag <- c(LFS = "grey95", Healthy = "#B2DF8A")
col_previous <- c(Healthy = "#B2DF8A", Yes = "#fb9a99", No = "grey95")
col_age <- c(Healthy = "#B2DF8A", Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_type <- c(Healthy = "#B2DF8A", LOF = "#FC8D62", "Splice Site"= "#66C2A5", "Missense 1" = "#E5C494", 
              "Missense 2" = "#FFD92F", "Missense 3" = "#8DA0CB", "Missense 5" = "#E78AC3")
col_metric <- c("Nucleosome Footprint" = "#E69F00", "Fragment Ratio" = "#0072B2", "Fragment Length" = "#D55E00", 
                "Accessibility TFBS" = "#CC79A7", "Dinucleotide" = "#6A3D9A", "Accessibility TCGA/DHS" = "#56B4E9")

### Set annotations
top_annotation <- HeatmapAnnotation(#"Metrics" = anno_barplot(t(data_class),
                                    #                         gp = gpar(fill = col_metric,
                                    #                                   col = col_metric),
                                    #                         bar_width = 1,
                                    #                         axis_param = list(gp = gpar(fontsize = 8)),
                                    #                         border = FALSE),
                                    "Integrated LFS\nFragmentation Score" = anno_barplot(data$integrated,
                                                                                        gp = gpar(fill = "black",
                                                                                                  col = "black"),
                                                                                        bar_width = 1,
                                                                                        axis_param = list(gp = gpar(fontsize = 8)),
                                                                                        border = FALSE),
                                    "Breakpoint" = anno_barplot(data$breakpoint,
                                                                        gp = gpar(fill = "grey65",
                                                                                  col = "grey65"),
                                                                        bar_width = 1,
                                                                        axis_param = list(gp = gpar(fontsize = 0)),
                                                                        border = FALSE),
                                    "Accessibility TFBS" = anno_barplot(data$TFBS,
                                                                        gp = gpar(fill = "#CC79A7",
                                                                                  col = "#CC79A7"),
                                                                        bar_width = 1,
                                                                        axis_param = list(gp = gpar(fontsize = 0)),
                                                                        border = FALSE),
                                    "Accessibility TCGA/DHS" = anno_barplot(data$TCGA_DHS,
                                                                            gp = gpar(fill = "#56B4E9",
                                                                                      col = "#56B4E9"),
                                                                            bar_width = 1,
                                                                            axis_param = list(gp = gpar(fontsize = 0)),
                                                                            border = FALSE),
                                    "Fragment Ratio" = anno_barplot(data$ratio,
                                                                    gp = gpar(fill = "#0072B2",
                                                                              col = "#0072B2"),
                                                                    bar_width = 1,
                                                                    axis_param = list(gp = gpar(fontsize = 0)),
                                                                    border = FALSE),
                                    "Nucleosome Footprint" = anno_barplot(data$peaks,
                                                                          gp = gpar(fill = "#E69F00",
                                                                                    col = "#E69F00"),
                                                                          bar_width = 1,
                                                                          axis_param = list(gp = gpar(fontsize = 0)),
                                                                          border = FALSE),
                                    "Fragment Length" = anno_barplot(data$insert_size,
                                                                     gp = gpar(fill = "#D55E00",
                                                                               col = "#D55E00"),
                                                                     bar_width = 1,
                                                                     axis_param = list(gp = gpar(fontsize = 0)),
                                                                     border = FALSE),
                                    "Dinucleotide" = anno_barplot(data$dinucleotide,
                                                                  gp = gpar(fill = "#CC79A7",
                                                                            col = "#CC79A7"),
                                                                  bar_width = 1,
                                                                  axis_param = list(gp = gpar(fontsize = 0)),
                                                                  border = FALSE),
                                    border = FALSE,
                                    annotation_height = unit(c(1.5, rep(0.5, 7)), "cm"),
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 10),
                                    gap = unit(2, "mm"),
                                    show_legend = FALSE)

bot_annotation <- HeatmapAnnotation("Germline Mutation" = data$mutation_type,
                                    "Cancer History" = data$previous_cancer,
                                    "Patient Age" = data$Age,
                                    "Sample Type" = data$ActualClass,
                                    col = list("Germline Mutation" = col_type,
                                               "Cancer History" = col_previous,
                                               "Patient Age" = col_age,
                                               "Sample Type" = col_diag),
                                    border = FALSE,
                                    height = unit(1.75, "cm"),
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 10),
                                    simple_anno_size = unit(0.3, "cm"),
                                    show_legend = FALSE)

### Make annotation legend
annotation_legend = packLegend(list = list(Legend(title = "Sample Type", 
                                                  at = c("Healthy Control", "TP53m-carrier"),
                                                  legend_gp = gpar(fill = c("#B2DF8A", "grey95")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Patient Age", 
                                                  at = c("Adult", "Pediatric"),
                                                  legend_gp = gpar(fill = c("#6A3D9A", "#CAB2D6")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer History",
                                                  at = c("Yes"),
                                                  legend_gp = gpar(fill = c("#fb9a99")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Germline Mutation",
                                                  at = c("LOF", "Splice Site", "Missense 1", "Missense 2", "Missense 3", "Missense 5"),
                                                  legend_gp = gpar(fill = c("#FC8D62", "#66C2A5", "#E5C494", "#FFD92F", "#8DA0CB", "#E78AC3")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 2),
                                           #Legend(title = "Metric",
                                          #        at = c("Nucleosome Footprint", "Fragment Ratio", "Fragment Length", 
                                          #               "Dinucleotide", "Accessibility TFBS", "Accessibility TCGA/DHS"),
                                          #        legend_gp = gpar(fill = c("#E69F00", "#0072B2", "#D55E00", "#6A3D9A", "#CC79A7", "#56B4E9")),
                                          #        title_gp = gpar(fontsize = 8),
                                          #        labels_gp = gpar(fontsize = 7),
                                          #        grid_height = unit(1, "mm")),
                                           Legend(title = "Score",
                                                  at = c(0, 0.5, 1),
                                                  col_fun = col_fun,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_width = unit(2, "mm"),
                                                  grid_height = unit(1, "mm"),
                                                  direction = "horizontal")),
                               direction = "horizontal")


### Set order
col_order <- colnames(data_class)

### Make heatmap
pdf(file.path(path, "integrated_scores_hbc_v_lfs.pdf"), width = 8, height = 5.5)
Heatmap <- Heatmap(data_class,
                   col = col_fun,
                   column_order = col_order,
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   top_annotation = top_annotation,
                   bottom_annotation = bot_annotation,
                   column_title_gp = gpar(fontsize = 6),
                   row_names_gp = gpar(fontsize = 10),
                   row_names_side = "left",
                   row_title_rot = 0,
                   row_dend_width = unit(0, "cm"),
                   border = FALSE,
                   height = unit(2, "cm"))
draw(Heatmap, show_annotation_legend = FALSE, heatmap_legend_list = annotation_legend, heatmap_legend_side = "bottom")
dev.off()

### Comparisons
data_comp <- data[data$ActualClass == "LFS", ]

data_stats_history <- data_comp %>%
  group_by(previous_cancer)%>% 
  dplyr::summarise(Median=median(integrated, na.rm = TRUE),
                   Mean=mean(integrated, na.rm = TRUE),
                   SD=sd(integrated, na.rm = TRUE),
                   N=n())
a <- t.test(data_comp$integrated[data_comp$previous_cancer=="No"], data_comp$integrated[data_comp$previous_cancer=="Yes"])$p.value
data_stats_history$pvalue <- c(1, a)
data_stats_history$annot <- ifelse(data_stats_history$pvalue < 0.05 & data_stats_history$pvalue > 0.01, "*",
                                  ifelse(data_stats_history$pvalue < 0.01 & data_stats_history$pvalue > 0.001, "**",
                                         ifelse(data_stats_history$pvalue < 0.001, "***", "")))

data_stats_age <- data_comp %>%
  group_by(Age)%>% 
  dplyr::summarise(Median=median(integrated, na.rm = TRUE),
                   Mean=mean(integrated, na.rm = TRUE),
                   SD=sd(integrated, na.rm = TRUE),
                   N=n())
a <- t.test(data_comp$integrated[data_comp$Age=="Pediatric"], data_comp$integrated[data_comp$Age=="Adult"])$p.value
data_stats_age$pvalue <- c(1, a)
data_stats_age$annot <- ifelse(data_stats_age$pvalue < 0.05 & data_stats_age$pvalue > 0.01, "*",
                                   ifelse(data_stats_age$pvalue < 0.01 & data_stats_age$pvalue > 0.001, "**",
                                          ifelse(data_stats_age$pvalue < 0.001, "***", "")))

data_stats_mut <- data_comp %>%
  group_by(mutation_type)%>% 
  dplyr::summarise(Median=median(integrated, na.rm = TRUE),
                   Mean=mean(integrated, na.rm = TRUE),
                   SD=sd(integrated, na.rm = TRUE),
                   N=n())
vars <- data_stats_mut$mutation_type
for (i in vars) {
  t_test <- c()
  for (var in vars) {
    a <- t.test(data_comp$integrated[data_comp$mutation_type==i], data_comp$integrated[data_comp$mutation_type==var])$p.value
    t_test <- c(t_test, a)
  }
  data_stats_mut <- cbind(data_stats_mut, t_test)
}
colnames(data_stats_mut)[6:11] <- as.vector(vars)

## Plot Figures
### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 12),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

data_comp <- data_comp %>%
  group_by(previous_cancer) %>%
  arrange(integrated)
fig_his <- ggplot(data_comp) +
  geom_boxplot(aes(previous_cancer, integrated, fill = previous_cancer), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(previous_cancer, integrated), position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_history, aes(x = previous_cancer, y = -0.05, label = N), size = 4) +
  geom_text(data = data_stats_history, aes(x = previous_cancer, y = 1.1, label = annot), size = 5) +
  ggtitle("Cancer\nHistory") + 
  xlab("") +
  ylab("Integrated Score") +
  scale_fill_manual(values = c("grey65", "red")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(-0.1, 1.2), expand = c(0,0))
fig_his

data_comp <- data_comp %>%
  group_by(Age) %>%
  arrange(integrated)
fig_age <- ggplot(data_comp) +
  geom_boxplot(aes(Age, integrated, fill = Age), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(Age, integrated), position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_age, aes(x = Age, y = -0.05, label = N), size = 4) +
  geom_text(data = data_stats_age, aes(x = Age, y = 1.1, label = annot), size = 5) +
  ggtitle("Age") + 
  xlab("") +
  ylab("Integrated Score") +
  scale_fill_manual(values = c("grey65", "red")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(-0.1, 1.2), expand = c(0,0))
fig_age

data_comp <- data_comp %>%
  group_by(mutation_type) %>%
  arrange(integrated)
fig_mut <- ggplot(data_comp, aes(mutation_type, integrated)) +
  geom_boxplot(aes(fill = mutation_type), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_mut, aes(x = mutation_type, y = -0.05, label = N), size = 4) +
  ggtitle("Germline Mutation Type") + 
  xlab("") +
  ylab("Integrated Score") +
  scale_fill_manual(values = c("#FB9A99", "#A6CEE3", "#FDBF6F","#CAB2D6", "grey65", "#B2DF8A")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(-0.1, 1.1), expand = c(0,0)) 

fig_mut

data_age <- data[data$ActualClass == "Healthy", ]
data_age <- merge(data_age[, c("sample", "integrated", "ActualClass")], samples_healthy, by.x = "sample", by.y = "sWGS")
data_age <- data_age[, c("sample", "integrated", "ActualClass", "age")]
data_comp <- data_comp[, c("sample", "integrated", "ActualClass", "years")]
colnames(data_comp) <- c("sample", "integrated", "ActualClass", "age")
data_age <- rbind(data_age, data_comp)
data_age$age <- as.numeric(data_age$age)
data_age <- data_age[complete.cases(data_age), ]

fig_cor <- ggplot(data_age, aes(age, integrated)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  stat_regline_equation(label.y = 0.15, label.x = 5, aes(label = ..rr.label..), size = 5) +
  ggtitle("") + 
  xlab("Age (years)") +
  ylab("Integrated Score") +
  facet_grid(.~ActualClass) +
  theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))
fig_cor

Figure <- ggarrange(fig_his, fig_age, fig_mut, fig_cor,  align = "hv", widths = c(1, 1, 2, 3.5), nrow = 1)
Figure
ggsave(file.path(path, "integrated_scores_comparison.pdf"), Figure, device = "pdf", width = 12, height = 4.5, units = "in")

### Further investigation
data_neg <- data[data$ActualClass == "LFS", ]
data_neg <- data_neg[!(data_neg$ext_ID %in% cancer_patients), ]

data_neg <- data_neg %>%
  group_by(Age, previous_cancer) %>%
  arrange(integrated)
fig_neg <- ggplot(data_neg) +
  geom_boxplot(aes(Age, integrated, fill = previous_cancer), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(Age, integrated, group = previous_cancer), position = position_dodge2(.5), size = 1.5, pch = 16) +
  ggtitle("Age") + 
  xlab("") +
  ylab("Integrated Score") +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(-0.1, 1.2), expand = c(0,0))
fig_neg

