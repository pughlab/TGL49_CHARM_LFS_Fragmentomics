library(tidyverse)
library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)

### Set paths
path <- "data/end_motifs"
outdir <- ""
healthy_path <- "hbc/end_motifs"

### Find paths
data <- list.files(path, "motifs.txt", full.names = TRUE)
data <- data[grepl("genome2", data)]
data_normal <- list.files(healthy_path, "motifs.txt", full.names = TRUE)
data_normal <- data_normal[grepl("genome2", data_normal)]

### Import data 
data <- read.delim(data)
data_normal <- read.delim(data_normal)
data_samples <- read.delim("sample_list.txt")

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
#abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]

### Calculate Healthy Median
hbc_median <- rowMedians(as.matrix(data_normal[,-1]))

### Combine dataframes for heatmap
healthy_samples <- as.data.frame(colnames(data_normal[, -1]))
colnames(healthy_samples) <- c("sWGS")
healthy_samples$Age <- "HBC"
healthy_samples$cancer_status <- "HBC"
healthy_samples$previous_cancer <- "HBC"
healthy_samples$mutation_type <- "HBC"

samples <- bind_rows(data_samples, healthy_samples)

data <- merge(data, data_normal, by = "motif")
row.names(data) <- data$motif
samples <- samples[samples$sWGS %in% colnames(data), ]
data <- data[, samples$sWGS]
data <- data/hbc_median
data <- as.matrix(data)

### Format metadata
samples$Age <- factor(samples$Age, levels = c("HBC", "pediatric", "adult"),
                      labels = c("Healthy", "Pediatric", "Adult"))
samples$cancer_status <- factor(samples$cancer_status, levels = c("HBC", "negative", "positive"),
                                labels = c("Healthy", "Negative", "Positive"))
samples$previous_cancer <- factor(samples$previous_cancer, levels = c("HBC", "yes", "no"),
                                  labels = c("Healthy", "Yes", "No"))
samples$mutation_type <- factor(samples$mutation_type,
                                levels = c("HBC", "LOF", "Splice", "1", "2", "3", "5"),
                                labels = c("Healthy", "LOF", "Splice", "1", "2", "3", "5"))

### Set clinical information
data_status <- as.matrix(samples$cancer_status)
row.names(data_status) <- samples$sWGS

data_age <- as.matrix(samples$Age)
row.names(data_age) <- samples$sWGS

data_type <- as.matrix(samples$mutation_type)
row.names(data_type) <- samples$sWGS

data_previous <- as.matrix(samples$previous_cancer)
row.names(data_previous) <- samples$sWGS

### Make Motif annotation
data_motif <- as.matrix(str_split_fixed(row.names(data), "", 4))
data_motif[data_motif == "A" | data_motif == "T"] <- "AT"
data_motif[data_motif == "G" | data_motif == "C"] <- "GC"
data_motif <- as.matrix(data_motif)
row.names(data_motif) <- row.names(data)

### Set annotation and heatmap colours
col_heat <- colorRamp2(c(0.5, 0.7, 0.9, 1.1, 1.3, 1.5), 
                       c("#fde725", "#7ad151", "#22a884", "#2a788e", "#414487", "#440154"))
col_status <- c(Healthy = "#B2DF8A", Positive = "#fb9a99", Negative = "#a6cee3")
col_previous <- c(Healthy = "#B2DF8A", Yes = "#fb9a99", No = "grey95")
col_age <- c(Healthy = "#B2DF8A", Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_type <- c(Healthy = "#B2DF8A", LOF = "#FC8D62", Splice = "#66C2A5", "1" = "#E5C494", "2" = "#FFD92F", "3" = "#8DA0CB", "5" = "#E78AC3")
col_motif <- c(AT = "#FC766AFF", GC = "#5B84B1FF")

### Set additional annotations
top_annotation <- HeatmapAnnotation("Patient Type" = data_age,
                                    "Cancer History" = data_previous,
                                    "Cancer Status" = data_status,
                                    "Mutation Type" = data_type,
                                    show_annotation_name = TRUE,
                                    border = FALSE,
                                    col = list("Cancer Status" = col_status, 
                                               "Patient Type" = col_age,
                                               "Cancer History" = col_previous,
                                               "Mutation Type" = col_type),
                                    annotation_name_gp = gpar(fontsize = 8),
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    simple_anno_size = unit(0.25, "cm"),
                                    show_legend = c(FALSE, FALSE, FALSE, FALSE))

right_annotation <- rowAnnotation("End Motif\n5' to 3'" = data_motif,
                                  border = FALSE,
                                  col = list("End Motif\n5' to 3'" = col_motif),
                                  annotation_name_gp = gpar(fontsize = 8),
                                  annotation_name_side = "top",
                                  annotation_name_rot = 0,
                                  simple_anno_size = unit(0.3, "cm"),
                                  show_legend = c(FALSE))

### Set legend labels
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
                                                  at = c("Healthy", "Negative", "Positive"),
                                                  legend_gp = gpar(fill = c("#B2DF8A", "#a6cee3", "#fb9a99")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Mutation Type",
                                                  at = c("LOF", "Splice", "Missense 1", "Missense 2", "Missense 3", "Missense 5"),
                                                  legend_gp = gpar(fill = c("#FC8D62", "#66C2A5", "#E5C494", "#FFD92F", "#8DA0CB", "#E78AC3")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "End Motif",
                                                  at = c("A/T", "G/C"),
                                                  legend_gp = gpar(fill = c("#FC766AFF", "#5B84B1FF")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Fold Change",
                                                  at = c(0.5, 1, 1.5),
                                                  col_fun = col_heat,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "vertical",
                                                  grid_width = unit(2, "mm"),
                                                  legend_height = unit(1.5, "cm"))),
                               direction = "horizontal")

## Generate Heatmap
pdf(file.path(outdir, "end_motifs_heatmap.pdf"), height = 6, width = 5)
Heatmap <- Heatmap(data,
                   col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   top_annotation = top_annotation,
                   right_annotation = right_annotation,
                   row_labels = NULL,
                   column_title_gp = gpar(fontsize = 6),
                   row_title = "5' Tetranucleotide Motifs",
                   row_title_gp = gpar(fontsize = 8),
                   row_title_rot = 90,
                   row_dend_width = unit(1.5, "cm"),
                   border = FALSE)
draw(Heatmap, heatmap_legend_list = annotation_legend, heatmap_legend_side = "bottom", show_annotation_legend = FALSE)
dev.off()

cluster <- draw(Heatmap)
cluster <- column_order(cluster)
cluster <- data[,cluster]
row.names(samples) <- samples$sWGS
cluster <- samples[colnames(cluster), ]



