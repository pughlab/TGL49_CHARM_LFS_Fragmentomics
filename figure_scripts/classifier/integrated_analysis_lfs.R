library(caret)
library(dplyr)
library(tidyverse)
library(pROC)
library(ggpubr)
library(ComplexHeatmap)
library(lemon)

### Set variables
path <- "data/classifier"

### Find files
filenames <- list.files(path, ".rds", full.names = TRUE)
filenames <- filenames[grepl("integrated_lfs.rds", filenames)]

### Read in file
data <- readRDS(filenames)
samples <- read.delim("sample_list.txt")
samples <- samples[!(samples$notes == "NS"), ]

ichor <- read.delim("data/CHARM_LFS_ichorCNA_summary_reviewed.txt")
hbc_scores <- read.delim(file.path(path, "longitudinal_hbc_scores.txt"))

### Unlist files
data <- unlist(data, recursive = FALSE)
data <- bind_rows(data)
data$ActualClass <- factor(data$ActualClass, levels = c("negative", "positive"))

### Aggregate scores per sample and calculate delta values
data_int <- merge(data, samples[, c("ext_ID", "timepoint", "Age", "cancer_status", "previous_cancer", "stage", "cancer_type", "sWGS")], by.x = "sample", by.y = "sWGS")
data_int <- data_int[data_int$ext_ID %in% samples$ext_ID[duplicated(samples$ext_ID)], ]

data_baselines <- data_int[data_int$timepoint == "0", ] %>%
  group_by(ext_ID) %>%
  dplyr::summarise(mean=mean(positive))

count <- data_int %>%
  group_by(ext_ID) %>%
  dplyr::summarise(N=n())

baselines <- rep(data_baselines$mean, count$N)

data_int <- data_int[order(data_int$ext_ID), ]
data_int$positive <- (data_int$positive - baselines)/baselines

data_int <- data_int %>%
  group_by(sample) %>%
  dplyr::summarise(score=mean(positive),
                   sd=sd(positive),
                   N=n())
data_int$CI <- data_int$sd/sqrt(data_int$N)*qt(0.975, data_int$N -1)

### Keep serial phenoconverters samples
data_int <- merge(data_int, samples[, c("ext_ID", "timepoint", "Age", "cancer_status", "previous_cancer", "stage", "cancer_type", "sWGS")], by.x = "sample", by.y = "sWGS")
data_int <- data_int[data_int$ext_ID %in% data_int$ext_ID[duplicated(data_int$ext_ID)], ]

### Merge with ichorCNA
ichor$ichor <- ifelse(ichor$TF > ichor$TF_short, ichor$TF, ichor$TF_short)
ichor$ichor <- as.numeric(ichor$ichor)
ichor <- ichor[, c("sample", "ichor")]
data_int <- merge(data_int, ichor, by = "sample")

### Calculate slopes
data_int$timepoint <- factor(data_int$timepoint, levels = c(0:10))
data_int <- data_int[order(data_int$ext_ID,
                           data_int$timepoint), ]
data_int$slope <- ave(data_int$score, data_int$ext_ID, FUN = function(x) c(0, diff(x)))

### Split phenoconverters and negatives
phenos <- data_int$ext_ID[data_int$cancer_status == "positive"]
phenos <- phenos[phenos %in% data_int$ext_ID[data_int$cancer_status == "negative"]]

data_neg <- data_int[data_int$cancer_status == "negative", ]
data_neg <- data_neg[data_neg$ext_ID %in% data_neg$ext_ID[duplicated(data_neg$ext_ID)], ]
data_neg <- data_neg[!(data_neg$ext_ID %in% phenos), ]

data_int <- data_int[data_int$ext_ID %in% phenos, ]
data_int <- data_int[order(data_int$ext_ID,
                           data_int$timepoint),]

### Make dataframes with timepoints for stats
hbc_scores <- hbc_scores[!(hbc_scores$timepoint == 0), ]

forwards <- data.frame(ext_ID = c("LFS45", 
                                  rep("LFS55", 2), 
                                  "LFS68", 
                                  "LFS71", 
                                  rep("LFS3", 6),
                                  rep("LFS5", 4),
                                  rep("LFS15", 3),
                                  rep("LFS34", 2),
                                  rep("LFS51", 2),
                                  rep("LFS59", 2),
                                  rep("LFS78"), 2),
                       timepoint = c(1, 
                                     1,2, 
                                     1, 
                                     3, 
                                     1,2,6,7,8,9,
                                     3,4,5,6,
                                     2,3,4,
                                     1,2,
                                     1,2,
                                     1,2,
                                     1,2))

revs <- data.frame(ext_ID = c(rep("LFS2", 5),
                              rep("LFS35", 2),
                              rep("LFS37", 3),
                              rep("LFS40", 3),
                              "LFS57", 
                              "LFS77", 
                              "LFS81", 
                              "LFS88", 
                              "LFS89"),
                   timepoint = c(1,2,3,4,5,
                                 1,2,
                                 1,2,3,
                                 1,2,3,
                                 2,
                                 1,
                                 1,
                                 1,
                                 1))

### Calculate delta statistics
score_neg <- data_neg[!(data_neg$timepoint == 0) & 
                        !(data_neg$ext_ID  %in% c("LFS28", "LFS45", "LFS55", "LFS68", "LFS71")), ]
stats_score_neg <- data.frame(mean = mean(score_neg$score),
                              sd = sd(score_neg$score))

score_forward <- merge(data_int, forwards, by = c("ext_ID", "timepoint"))
stats_score_forward <- data.frame(mean = mean(score_forward$score),
                                  sd = sd(score_forward$score))

score_rev <- merge(data_int, revs, by = c("ext_ID", "timepoint"))
stats_score_reverse <- data.frame(mean = mean(score_rev$score),
                                  sd = sd(score_rev$score))

stats_score_p <- data.frame(comp = c("Neg", "Forward", "Reverse"),
                            pvalue = c(t.test(hbc_scores$score, score_neg$score)$p.value,
                                       t.test(hbc_scores$score, score_forward$score)$p.value,
                                       t.test(hbc_scores$score, score_rev$score)$p.value))

### Calculate slope statistics
slope_neg <- data_neg[!(data_neg$timepoint == 0) & 
                        !(data_neg$ext_ID  %in% c("LFS28", "LFS45", "LFS55", "LFS68", "LFS71")), ]
stats_slope_neg <- data.frame(mean = mean(slope_neg$slope),
                              sd = sd(slope_neg$slope))

slope_forward <- merge(data_int, forwards, by = c("ext_ID", "timepoint"))
stats_slope_forward <- data.frame(mean = mean(slope_forward$slope),
                                  sd = sd(slope_forward$slope))

slope_rev <- merge(data_int, revs, by = c("ext_ID", "timepoint"))
stats_slope_reverse <- data.frame(mean = mean(slope_rev$slope),
                                  sd = sd(slope_rev$slope))

stats_slope_p <- data.frame(comp = c("Neg", "Forward", "Reverse"),
                            pvalue = c(t.test(hbc_scores$slope, score_neg$slope)$p.value,
                                       t.test(hbc_scores$slope, score_forward$slope)$p.value,
                                       t.test(hbc_scores$slope, score_rev$slope)$p.value))

### Set samples to graph and order factors (phenoconverters)
data_int$cancer_status <- factor(data_int$cancer_status,
                                levels = c("negative", "positive"),
                                labels = c("Negative", "Positive"))
data_int$Age <- factor(data_int$Age,
                      levels = c("adult", "pediatric"),
                      labels = c("Adult", "Pediatric"))
data_int$previous_cancer <- factor(data_int$previous_cancer, levels = c("yes", "no", ""),
                                  labels = c("Yes", "No", "Unknown"))
data_int$stage <- factor(data_int$stage, levels = c("high", "low", "", "unknown"),
                            labels = c("Stage III/IV", "Stage 0/I/II", "none", "Unknown"))
data_int$cancer_type <- factor(data_int$cancer_type, levels = c("adrenal", "appendiceal_adenocarcinoma", "astrocytoma", "bladder", "breast", "endometrial", 
                                                                "glioma", "leukemia", "lung", "lymphoma", "osteosarcoma", "prostate", "sarcoma", ""),
                               labels = c("Adrenal", "Adenocarcinoma", "Glioma", "Bladder", "Breast", "Endometrial", "Glioma",
                                          "Leukemia", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma", "NA"))
data_matrix <- matrix(data = 0, ncol = nrow(data_int))
data_int$ext_ID <- gsub("LFS", "", data_int$ext_ID)
colnames(data_matrix) <- data_int$sample

### Set samples to graph and order factors (negative serial)
data_neg$cancer_status <- factor(data_neg$cancer_status,
                                 levels = c("negative", "positive"),
                                 labels = c("Negative", "Positive"))
data_neg$Age <- factor(data_neg$Age,
                       levels = c("adult", "pediatric"),
                       labels = c("Adult", "Pediatric"))
data_neg$previous_cancer <- factor(data_neg$previous_cancer, levels = c("yes", "no", ""),
                                   labels = c("Yes", "No", "Unknown"))
data_neg$developed <- ifelse(data_neg$ext_ID %in% c("LFS28", "LFS45", "LFS55", "LFS68", "LFS71"), "yes", "no")

data_neg <- data_neg[order(data_neg$developed,
                           factor(data_neg$ext_ID, levels = paste0("LFS", 1:71)),
                           data_neg$timepoint),]
data_matrix2 <- matrix(data = 0, ncol = nrow(data_neg))
data_neg$ext_ID <- gsub("LFS", "", data_neg$ext_ID)
colnames(data_matrix2) <- data_neg$sample

### Set annotations (phenoconverters)
# Make fragmentation score
data_score <- data.frame(score = data_int$score,
                         min = data_int$score - data_int$sd,
                         max = data_int$score + data_int$sd,
                         zero = 0,
                         upper = 0.57)
data_score <- as.matrix(data_score)
row.names(data_score) <- data_int$sample

# Slope
data_slope <- data.frame(slope = data_int$slope,
                         zero = 0,
                         upper = 0.569)
data_slope <- as.matrix(data_slope)
row.names(data_slope) <- data_int$sample

# Set Age
data_age <- as.matrix(data_int$Age)
row.names(data_age) <- data_int$sample

# Set Cancer status
data_cancer <- as.matrix(data_int$cancer_status)
row.names(data_cancer) <- data_int$sample

# Set Previous cancers
data_history <- as.matrix(data_int$previous_cancer)
row.names(data_history) <- data_int$sample

# Set Cancer Stage
data_stage <- as.matrix(data_int$stage)
row.names(data_stage) <- data_int$sample

# Set Cancer Type
data_type <- as.matrix(data_int$cancer_type)
row.names(data_type) <- data_int$sample

# ichorCNA
data_ichor <- as.matrix(data_int$ichor)
data_ichor <- cbind(data_ichor, 0.03)
row.names(data_ichor) <- data_int$sample

### Set annotations (serial negatives)
# Make fragmentation score
data_score2 <- data.frame(score = data_neg$score,
                          min = data_neg$score - data_neg$sd,
                          max = data_neg$score + data_neg$sd,
                          zero = 0,
                          upper = 0.57)
data_score2 <- as.matrix(data_score2)
row.names(data_score2) <- data_neg$sample

# Slope
data_slope2 <- data.frame(slope = data_neg$slope,
                          zero = 0,
                          upper = 0.569)
data_slope2 <- as.matrix(data_slope2)
row.names(data_slope2) <- data_neg$sample

# Set Age
data_age2 <- as.matrix(data_neg$Age)
row.names(data_age2) <- data_neg$sample

# Set cancer development
data_dev <- as.matrix(data_neg$developed)
row.names(data_dev) <- data_neg$sample

# Set Previous cancers
data_history2 <- as.matrix(data_neg$previous_cancer)
row.names(data_history2) <- data_neg$sample

# Set ichor
data_ichor2 <- as.matrix(data_neg$ichor)
data_ichor2 <- cbind(data_ichor2, 0.03)
row.names(data_ichor2) <- data_neg$sample

### Set colours
col_cancer <- c(Positive = "#fb9a99", Negative = "#a6cee3")
col_age <- c(Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_history <- c(Yes = "#fb9a99", No = "grey95", Unknown = "grey95")
col_dev <- c(yes = "#fb9a99", no = "grey95")
col_stage <- c("Stage III/IV" = "#fb9a99", "Stage 0/I/II" = "#ebdc78", "none" = "grey95", Unknown = "grey65")
col_type <- c(Adenocarcinoma = "#fd7f6f", Adrenal = "#7eb0d5", Bladder = "#b2e061", Breast = "#bd7ebe", Endometrial = "#ffb55a", 
              Leukemia = "#add3ac", Glioma = "grey65", Lung = "#fdcce5", Lymphoma = "#8bd3c7", Osteosarcoma = "#ebdc78", 
              Prostate = "#b3d4ff", Sarcoma = "#beb9db", Unknown = "grey95", "NA" = "grey95")

## Set annotations
bot_annotation <- HeatmapAnnotation("Cancer Type" = data_type,
                                    "Cancer Stage" = data_stage,
                                    "Cancer History" = data_history,
                                    "Patient Type" = data_age,
                                    border = FALSE,
                                    annotation_name_side = "left",
                                    show_annotation_name = TRUE,
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 8),
                                    col = list("Patient Type" = col_age,
                                               "Cancer History" = col_history,
                                               "Cancer Status" = col_cancer,
                                               "Cancer Type" = col_type,
                                               "Cancer Stage" = col_stage),
                                    show_legend = FALSE,
                                    simple_anno_size = unit(0.3, "cm"))

top_annotation <- HeatmapAnnotation("Cancer Status" = data_cancer,
                                    "Cancer Fragmentation\nDelta (CFD)" = anno_lines(data_score,
                                                       ylim = c(-1, 2.5),
                                                       add_points = TRUE,
                                                       pch = 16,
                                                       gp = gpar(col = c("black", "red", "red", "black", "red"),
                                                                 lty = c("solid", "dashed", "dashed", "dashed", "dashed"),
                                                                 alpha = c(1, 0.25, 0.25, 0.5, 1)),
                                                       pt_gp = gpar(col = c("black")),
                                                       size = unit(c(1,0,0,0,0), "mm"),
                                                       axis_param = list(side = "left",
                                                                         labels_rot = 0,
                                                                         gp = gpar(fontsize = 8)),
                                                       border = FALSE),
                                    "Cancer Fragmentation\nSlope (CFS)" = anno_lines(data_slope,
                                                       ylim = c(-1, 1.25),
                                                       add_points = TRUE,
                                                       pch = 16,
                                                       gp = gpar(col = c("black", "black", "red"),
                                                                 lty = c("solid", "dashed", "dashed"),
                                                                 alpha = c(1, 0.5, 1)),
                                                       pt_gp = gpar(col = c("black")),
                                                       size = unit(c(1,0,0), "mm"),
                                                       axis_param = list(side = "left",
                                                                         labels_rot = 0,
                                                                         gp = gpar(fontsize = 8)),
                                                       border = FALSE),
                                    "Tumor Fraction\n(Copy Number)" = anno_lines(data_ichor,
                                                                   add_points = TRUE,
                                                                   pch = 16,
                                                                   gp = gpar(col = c("black", "black"),
                                                                             lty = c("solid", "dashed"),
                                                                             alpha = c(1, 0.5)),
                                                                   size = unit(c(1,0), "mm"),
                                                                   axis_param = list(side = "left",
                                                                                     labels_rot = 0,
                                                                                     gp = gpar(fontsize = 8)),
                                                                   border = FALSE),
                                    gap = unit(2, "mm"),
                                    border = FALSE,
                                    col = list("Cancer Status" = col_cancer),
                                    show_annotation_name = TRUE,
                                    annotation_name_rot = 0,
                                    annotation_name_side = "left",
                                    annotation_name_gp = gpar(fontsize = 8),
                                    show_legend = FALSE,
                                    height = unit(5.3, "cm"),
                                    simple_anno_size = unit(0.3, "cm"))

bot_annotation2 <- HeatmapAnnotation("Developed Cancer" = data_dev,
                                     "Cancer History" = data_history2,
                                     "Patient Type" = data_age2,
                                     border = FALSE,
                                     annotation_name_side = "left",
                                     show_annotation_name = TRUE,
                                     annotation_name_rot = 0,
                                     annotation_name_gp = gpar(fontsize = 8),
                                     col = list("Developed Cancer" = col_dev,
                                                "Patient Type" = col_age,
                                                "Cancer History" = col_history),
                                     show_legend = FALSE,
                                     simple_anno_size = unit(0.3, "cm"))

top_annotation2 <- HeatmapAnnotation("CFD" = anno_lines(data_score2,
                                                        ylim = c(-1, 2.5),
                                                        add_points = TRUE,
                                                        pch = 16,
                                                        gp = gpar(col = c("black", "red", "red", "black", "red"),
                                                                  lty = c("solid", "dashed", "dashed", "dashed", "dashed"),
                                                                  alpha = c(1, 0.25, 0.25, 0.5, 1)),
                                                        pt_gp = gpar(col = c("black")),
                                                        size = unit(c(1,0,0,0,0), "mm"),
                                                        axis_param = list(side = "left",
                                                                          labels_rot = 0,
                                                                          gp = gpar(fontsize = 8)),
                                                        border = FALSE),
                                     "CFS" = anno_lines(data_slope2,
                                                        ylim = c(-1, 1.25),
                                                        add_points = TRUE,
                                                        pch = 16,
                                                        gp = gpar(col = c("black", "black", "red"),
                                                                  lty = c("solid", "dashed", "dashed"),
                                                                  alpha = c(1, 0.5, 1)),
                                                        pt_gp = gpar(col = c("black")),
                                                        size = unit(c(1,0,0), "mm"),
                                                        axis_param = list(side = "left",
                                                                          labels_rot = 0,
                                                                          gp = gpar(fontsize = 8)),
                                                        border = FALSE),
                                     "Tumor Fraction\n(Copy Number)" = anno_lines(data_ichor2,
                                                                                  ylim = c(0, 0.2),
                                                                                  add_points = TRUE,
                                                                                  pch = 16,
                                                                                  gp = gpar(col = c("black", "black"),
                                                                                            lty = c("solid", "dashed"),
                                                                                            alpha = c(1, 0.5)),
                                                                                  size = unit(c(1,0), "mm"),
                                                                                  axis_param = list(side = "left",
                                                                                                    labels_rot = 0,
                                                                                                    gp = gpar(fontsize = 8)),
                                                                                  border = FALSE),
                                     gap = unit(2, "mm"),
                                     border = FALSE,
                                     show_annotation_name = TRUE,
                                     annotation_name_rot = 0,
                                     annotation_name_side = "left",
                                     annotation_name_gp = gpar(fontsize = 8),
                                     show_legend = FALSE,
                                     height = unit(5, "cm"))

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Patient Type", 
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
                                           Legend(title = "Cancer Status",
                                                  at = c("Positive", "Negative"),
                                                  legend_gp = gpar(fill = c("#fb9a99", "#a6cee3")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Stage",
                                                  at = c("Stage 0/I/II", "Stage III/IV", "Unknown/None"),
                                                  legend_gp = gpar(fill = c("#ebdc78", "#fb9a99", "grey65")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Cancer Type",
                                                  at = c("Adenocarcinoma", "Adrenal", "Bladder", "Breast", "Endometrial", "Leukemia", 
                                                         "Glioma", "Lung", "Lymphoma", "Osteosarcoma", "Prostate", "Sarcoma"),
                                                  legend_gp = gpar(fill = c("#fd7f6f", "#7eb0d5", "#b2e061","#bd7ebe","#ffb55a", "#add3ac",
                                                                            "grey65", "#fdcce5","#8bd3c7","#ebdc78","#b3d4ff","#beb9db")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 3),
                                           Legend(title = "Metrics",
                                                  type = "lines",
                                                  at = c("Score", "95% CI", "Detection Limit"),
                                                  legend_gp = gpar(col = c("black", "red", "red"),
                                                                   alpha = c(1, 0.25, 1),
                                                                   lty = c("solid", "dashed", "dashed")),
                                                  background = "white",
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"))),
                               direction = "horizontal")

annotation_legend2 = packLegend(list = list(Legend(title = "Patient Type", 
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
                                            Legend(title = "Developed Cancer",
                                                   at = c("Yes"),
                                                   legend_gp = gpar(fill = c("#fb9a99")),
                                                   title_gp = gpar(fontsize = 8),
                                                   labels_gp = gpar(fontsize = 7),
                                                   grid_height = unit(1, "mm")),
                                            Legend(title = "Metrics",
                                                   type = "lines",
                                                   at = c("Score", "95% CI", "Detection Limit"),
                                                   legend_gp = gpar(col = c("black", "red", "red"),
                                                                    alpha = c(1, 0.25, 1),
                                                                    lty = c("solid", "dashed", "dashed")),
                                                   background = "white",
                                                   title_gp = gpar(fontsize = 8),
                                                   labels_gp = gpar(fontsize = 7),
                                                   grid_height = unit(1, "mm"))),
                                direction = "horizontal")

## Set orders
col_order <- colnames(data_matrix)
col_split <- data_int$ext_ID
order <- c("3","5", "15","34","51","59","78","2","35","37","40","57","77","81","88","89")
col_split <- factor(col_split, levels = order)

col_order2 <- colnames(data_matrix2)
col_split2 <- data_neg$ext_ID
col_split2 <- factor(col_split2, levels = unique(data_neg$ext_ID))

## Generate heatmap (phenoconverters)
pdf(file.path(path, "integration_cancer_status_pheno_delta3.pdf"), width = 8, height = 4)
oncoPrint <- Heatmap(data_matrix,
                     show_heatmap_legend = FALSE,
                     column_order = col_order,
                     top_annotation = top_annotation,
                     bottom_annotation = bot_annotation,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     column_split = col_split,
                     column_gap = unit(0.35, "mm"),
                     column_title_rot = 0,
                     column_title_gp = gpar(fontsize = 7),
                     column_names_side = "bottom",
                     column_names_rot = 90,
                     column_names_gp = gpar(fontsize = 7),
                     border = FALSE,
                     border_gp = gpar(col = "black"),
                     height = unit(0, "cm"),
                     width = unit(0.25*ncol(data_matrix), "cm"))
draw(oncoPrint, show_annotation_legend = FALSE, heatmap_legend_list = annotation_legend, heatmap_legend_side = "bottom")
dev.off()

## Generate heatmap (serial negs)
pdf(file.path(path, "integration_cancer_status_negs_delta3.pdf"), width = 6, height = 4)
oncoPrint <- Heatmap(data_matrix2,
                     show_heatmap_legend = FALSE,
                     column_order = col_order2,
                     top_annotation = top_annotation2,
                     bottom_annotation = bot_annotation2,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     column_split = col_split2,
                     column_gap = unit(0.35, "mm"),
                     column_title_rot = 0,
                     column_title_gp = gpar(fontsize = 7),
                     column_names_side = "bottom",
                     column_names_rot = 90,
                     column_names_gp = gpar(fontsize = 7),
                     border = FALSE,
                     border_gp = gpar(col = "black"),
                     height = unit(0, "cm"),
                     width = unit(0.25*ncol(data_matrix2), "cm"))
draw(oncoPrint, show_annotation_legend = FALSE, heatmap_legend_list = annotation_legend2, heatmap_legend_side = "bottom")
dev.off()

write.table(data_neg, file.path(path, "integration_negative_scores3.txt"), sep = "\t", row.names = FALSE)

### Separate multiple cancers to build timeliness
data_mult <- data_int[data_int$ext_ID %in% c("3", "5", "40"), ]
data_mult <- merge(data_mult, samples[, c("date_blood", "sWGS")], by.x = "sample", by.y = "sWGS")

### Convert dates to days
data_mult$date_blood <- as.Date(data_mult$date_blood)
data_date <- data_mult$date_blood[data_mult$timepoint == "0"]
counts <- data_mult %>%
  group_by(ext_ID) %>%
  dplyr::summarise(N=n())
data_date <- rep(data_date, c(10,8,5))
data_date <- as.Date(data_date)

data_mult$days <- difftime(data_mult$date_blood, data_date, units = "days")
data_mult$days <- as.numeric(gsub(" days", "", data_mult$days))/(365/12)

### Format for plotting
data_mult$ext_ID <- paste0("LFS", data_mult$ext_ID)
data_mult$ext_ID <- factor(data_mult$ext_ID, levels = c("LFS3", "LFS5", "LFS40"))

theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.position = "bottom",
               legend.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(size = 12),
               axis.text = element_text(size = 12),
               axis.title = element_text(size = 12),
               axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

plot <- ggplot(data_mult, aes(days, score)) +
  geom_smooth(se = FALSE, alpha = 0.5, color = "grey65") +
  geom_point(aes(color = cancer_status)) +
  geom_errorbar(aes(ymin = score - CI, ymax = score + CI, color = cancer_status), width = 0, size = 1) +
  scale_color_manual(values = c(Negative = "black", Positive = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = -1, size = 1) +
  xlab("Months from Baseline") +
  ylab("Integrated Score (Delta)") +
  labs(color = "Cancer Status") +
  facet_grid(.~ext_ID, scales = "free_x", space = "free_x") +
  facet_rep_grid(.~ext_ID, scales = "free_x", space = "free_x") +
  theme +
  scale_y_continuous(limits = c(-1.5, 1.25), breaks = c(-0.5,0,0.5,1))
plot

ggsave(file.path(path, "multicancer_patients3.pdf"), plot, width = 10, height = 4)



