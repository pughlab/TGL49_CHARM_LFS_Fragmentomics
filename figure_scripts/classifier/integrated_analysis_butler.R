library(dplyr)
library(tidyverse)
library(ggpubr)
library(lemon)

### Set variables
path <- "data/classifier"

### Find files
filenames <- list.files(path, ".rds", full.names = TRUE)
filenames <- filenames[grepl("integrated", filenames)]
filenames <- filenames[grepl("butler", filenames)]

type <- "lfs"
analysis <- "Cancer Fragmentation Score"
i <- 1

### Read in file
file <- filenames[[i]]

data <- readRDS(file)
samples <- read.delim("butler/sample_list_butler.txt")
samples <- samples[!(samples$patient == "HD-130"), ]

### Unlist files
data <- unlist(data, recursive = FALSE)
data <- bind_rows(data)
data$sample <- gsub("\\.", "-", data$sample)

### Aggregate scores per sample and calculate delta values
data_int <- merge(data, samples, by.x = "sample", by.y = "sWGS")

data_baselines <- data_int[data_int$timepoint == "0", ] %>%
  group_by(patient) %>%
  dplyr::summarise(mean=mean(positive))

count <- data_int %>%
  group_by(patient) %>%
  dplyr::summarise(N=n())

baselines <- rep(data_baselines$mean, count$N)

data_int$positive <- (data_int$positive - baselines)/baselines

data_int <- data_int %>%
  group_by(sample) %>%
  dplyr::summarise(score=mean(positive),
                   sd=sd(positive),
                   N=n())
data_int$CI <- data_int$sd/sqrt(data_int$N)*qt(0.975, data_int$N - 1)

### Format plotting data
data_int <- merge(data_int, samples, by.x = "sample", by.y = "sWGS")
data_int$patient <- gsub("HD-", "", data_int$patient)

### Calculate standard deviation
sd <- sd(data_int$score)
median <- median(data_int$score[!(data_int$timepoint==0)])
upper <- quantile(data_int$score, 0.99)
stats_score <- as.data.frame(cbind(median, sd, upper))

### Calculate slopes
data_int$timepoint <- factor(data_int$timepoint, levels = c(0:5))
data_int <- data_int[order(data_int$patient,
                           data_int$timepoint), ]
data_int$slope <- ave(data_int$score, data_int$patient, FUN = function(x) c(0, diff(x)))

data_slope <- data_int[!(data_int$timepoint == 0), ]
stats_slope <- data.frame(median = median(data_slope$slope),
                          sd = sd(data_slope$slope))

slope_upper <- quantile(data_slope$slope, 0.99)

### Make plotting table
data_matrix <- matrix(data = 0, ncol = nrow(data_int))

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

### Set annotations
top_annotation <- HeatmapAnnotation("Cancer Fragmentation\nDelta (CFD)" = anno_lines(data_score,
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
                                    gap = unit(2, "mm"),
                                    border = FALSE,
                                    show_annotation_name = TRUE,
                                    annotation_name_rot = 0,
                                    annotation_name_side = "left",
                                    annotation_name_gp = gpar(fontsize = 8),
                                    show_legend = FALSE,
                                    height = unit(3, "cm"))

## Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Metrics",
                                                  type = "lines",
                                                  at = c("Score", "95% CI", "Detection Limit"),
                                                  legend_gp = gpar(col = c("black", "red", "red"),
                                                                   alpha = c(1, 0.25, 1),
                                                                   lty = c("solid", "dashed", "dashed")),
                                                  background = "white",
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"))))

## Set orders
col_order <- colnames(data_matrix)
col_split <- data_int$patient
order <- unique(col_split)
col_split <- factor(col_split, levels = order)

### Make plot
pdf(file.path(path, "longitudinal_hbc.pdf"), width = 5.75, height = 1.75)
oncoPrint <- Heatmap(data_matrix,
                     show_heatmap_legend = FALSE,
                     column_order = col_order,
                     top_annotation = top_annotation,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     column_split = col_split,
                     column_gap = unit(0.35, "mm"),
                     column_title_rot = 0,
                     column_title_gp = gpar(fontsize = 7),
                     column_names_side = "bottom",
                     column_names_rot = 90,
                     column_names_gp = gpar(fontsize = 7),
                     border = TRUE,
                     border_gp = gpar(col = "black"),
                     height = unit(0, "cm"),
                     width = unit(0.25*ncol(data_matrix), "cm"))
draw(oncoPrint, show_annotation_legend = FALSE, heatmap_legend_list = annotation_legend, heatmap_legend_side = "right")
dev.off()

write.table(data_int, file.path(path, "longitudinal_hbc_scores.txt"), sep = "\t", row.names = FALSE)
