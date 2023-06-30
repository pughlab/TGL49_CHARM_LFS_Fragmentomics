library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)

### Set paths
path <- "data/griffin/TFBS"
outdir <- ""
healthy_path <- "hbc/griffin/TFBS"

data_griffin <- list.files(path, "features", full.names = TRUE)
data_samples <- "sample_list.txt"
data_normal <- list.files(healthy_path, "features", full.names = TRUE)

### Import data 
datalist <- lapply(data_griffin, function(x){read.delim(file = x)})
data_griffin <- do.call(rbind, datalist)
datalist <- lapply(data_normal, function(x){read.delim(file = x)})
data_normal <- do.call(rbind, datalist)
data_samples <- read.delim(data_samples)

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG", "TGL49_0019_Cf_U_PE_317_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude, abnormal)), ]
data_samples <- data_samples[!(data_samples$germline_mutation == "unknown"), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_griffin), ]
#data_samples <- data_samples[data_samples$sWGS %in% coverage$sWGS, ]

### Format griffin data
row.names(data_griffin) <- data_griffin$features
data_griffin <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS]

row.names(data_normal) <- data_normal$features
data_normal <- data_normal[, -1]

### Center midpoints and remove means
## Griffin
data_mean <- data_griffin[row.names(data_griffin) %like% "coverage", ]
data_mid <- data_griffin[row.names(data_griffin) %like% "mid", ]
data_amp <- data_griffin[row.names(data_griffin) %like% "amp", ]

data_mid <- data_mid - data_mean + 1
data_griffin <- bind_rows(data_mid, data_amp)

## Normals
data_mean <- data_normal[row.names(data_normal) %like% "coverage", ]
data_mid <- data_normal[row.names(data_normal) %like% "mid", ]
data_amp <- data_normal[row.names(data_normal) %like% "amp", ]

data_mid <- data_mid - data_mean + 1
data_normal <- bind_rows(data_mid, data_amp)

### Make healthy median
normal_median <- rowMedians(as.matrix(data_normal, na.rm = TRUE))
normal_sd <- rowSds(as.matrix(data_normal, na.rm = TRUE))

### Make matrix for healthy median
normal_medians <- as.matrix(normal_median)
row.names(normal_medians) <- row.names(data_normal)
normal_sds <- as.matrix(normal_sd)
row.names(normal_sds) <- row.names(data_normal)

### Perform differential analysis (Format samples and data)
data_de <- bind_cols(data_griffin, data_normal)
data_de <- as.data.frame(t(data_de))

samples_de <- as.data.frame(cbind(data_samples[, c("sWGS", "mutation_type", "cancer_status")], "LFS"))
colnames(samples_de) <- c("sWGS", "type", "cancer_status", "diag")
normal_de <- as.data.frame(cbind(colnames(data_normal), "HBC", "HBC", "HBC"))
colnames(normal_de) <- c("sWGS", "type", "cancer_status", "diag")
samples_de <- bind_rows(samples_de, normal_de)

data_de <- merge(data_de, samples_de, by.x = "row.names", by.y = "sWGS")
data_de_melt <- reshape2::melt(data_de, id = c("Row.names", "diag", "cancer_status", "type"))

### Calculate t-tests
t_test <- c()
vars <- unique(data_de_melt$variable)
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$diag == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$diag == "LFS" & data_de_melt$variable == var])$p.value
  t_test <- c(t_test, a)
}
t_test <- p.adjust(t_test, method = "BH")

t_test2 <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "2" & data_de_melt$variable == var])$p.value
  t_test2 <- c(t_test2, a)
}
t_test2 <- p.adjust(t_test2, method = "BH")

t_test3 <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "3" & data_de_melt$variable == var])$p.value
  t_test3 <- c(t_test3, a)
}
t_test3 <- p.adjust(t_test3, method = "BH")

t_testlof <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "LOF" & data_de_melt$variable == var])$p.value
  t_testlof <- c(t_testlof, a)
}
t_testlof <- p.adjust(t_testlof, method = "BH")

t_testsp <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "Splice" & data_de_melt$variable == var])$p.value
  t_testsp <- c(t_testsp, a)
}
t_testsp <- p.adjust(t_testsp, method = "BH")

t_test_neg <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$cancer_status == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$cancer_status == "negative" & data_de_melt$variable == var])$p.value
  t_test_neg <- c(t_test_neg, a)
}
t_test_neg <- p.adjust(t_test_neg, method = "BH")

t_test_pos <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$cancer_status == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$cancer_status == "positive" & data_de_melt$variable == var])$p.value
  t_test_pos <- c(t_test_pos, a)
}
t_test_pos <- p.adjust(t_test_pos, method = "BH")

t_test_lfs <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$cancer_status == "negative" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$cancer_status == "positive" & data_de_melt$variable == var])$p.value
  t_test_lfs <- c(t_test_lfs, a)
}
t_test_lfs <- p.adjust(t_test_lfs, method = "BH")

t_test <- as.data.frame(cbind(as.character(vars), t_test, t_test2, t_test3, t_testlof, t_testsp, t_test_neg, t_test_pos, t_test_lfs))
t_test <- as.data.frame(cbind(as.character(vars), sapply(t_test[, 2:ncol(t_test)], as.numeric)))

### Calculate change between negative vs positive patients
data_neg <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$cancer_status == "negative"]]
data_neg <- (data_neg - normal_median)/normal_sd
data_neg <- rowMedians(as.matrix(data_neg), na.rm = TRUE)
data_neg <- as.matrix(data_neg)
row.names(data_neg) <- row.names(data_griffin)

data_pos <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$cancer_status == "positive"]]
data_pos <- (data_pos - normal_median)/normal_sd
data_pos <- rowMedians(as.matrix(data_pos), na.rm = TRUE)
data_pos <- as.matrix(data_pos)
row.names(data_pos) <- row.names(data_griffin)

median_negative <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$cancer_status == "negative"]]))
sd_negative <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$cancer_status == "negative"]]))

data_lfs <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$cancer_status == "positive"]]
data_lfs <- (data_lfs - median_negative)/sd_negative
data_lfs <- rowMedians(as.matrix(data_lfs), na.rm = TRUE)
data_lfs <- as.matrix(data_lfs)
row.names(data_lfs) <- row.names(data_griffin)

### Calculate the distance from the healthy median
data_change <- (data_griffin - normal_median)/normal_sd
data_change2 <- (data_normal - normal_median)/normal_sd

### Set Median change of LFS cohort from HBCs and within mutation groups
data_freq <- rowMedians(as.matrix(data_change), na.rm = TRUE)
data_freq <- as.matrix(data_freq)
row.names(data_freq) <- row.names(data_change)

data_c2 <- data_change[, colnames(data_change) %in% data_samples$sWGS[data_samples$mutation_type == 2]]
data_freq2 <- rowMedians(as.matrix(data_c2), na.rm = TRUE)
data_freq2 <- as.matrix(data_freq2)
row.names(data_freq2) <- row.names(data_change)

data_c3 <- data_change[, colnames(data_change) %in% data_samples$sWGS[data_samples$mutation_type == 3]]
data_freq3 <- rowMedians(as.matrix(data_c3), na.rm = TRUE)
data_freq3 <- as.matrix(data_freq3)
row.names(data_freq3) <- row.names(data_change)

data_lof <- data_change[, colnames(data_change) %in% data_samples$sWGS[data_samples$mutation_type == "LOF"]]
data_freqlof <- rowMedians(as.matrix(data_lof), na.rm = TRUE)
data_freqlof <- as.matrix(data_freqlof)
row.names(data_freqlof) <- row.names(data_change)

data_sp <- data_change[, colnames(data_change) %in% data_samples$sWGS[data_samples$mutation_type == "Splice"]]
data_freqsp <- rowMedians(as.matrix(data_sp), na.rm = TRUE)
data_freqsp <- as.matrix(data_freqsp)
row.names(data_freqsp) <- row.names(data_change)

data_freq <- bind_cols(data_freq, data_freq2, data_freq3, data_freqlof, data_freqsp, data_neg, data_pos, data_lfs)
data_freq <- as.data.frame(data_freq)
row.names(data_freq) <- row.names(data_change)
colnames(data_freq) <- c("all", "c2", "c3", "lof", "splice", "negative", "positive", "cancer_status")

### Make Differential table and plot
de_table <- merge(data_freq, t_test, by.x = "row.names", by.y = "V1")
de_table <- merge(de_table, normal_medians, by.x = "Row.names", by.y = "row.names")
de_table <- merge(de_table, normal_sds, by.x = "Row.names", by.y = "row.names")
colnames(de_table) <- c("feature", "cohort", "c2", "c3", "lof", "splice", "negative", "positive", "cancer_status", "padj_cohort", "padj_c2", 
                        "padj_c3", "padj_lof", "padj_splice", "padj_negative", "padj_positive", "padj_cancer_status", "healthy_median", "healthy_sd")
write.table(de_table, file.path(outdir, "fragment_nucleosome_TFBS_differential.txt"), sep = "\t", row.names = FALSE)

### combine data
data <- bind_cols(data_change, data_change2)
samples <- as.data.frame(colnames(data_normal))
colnames(samples) <- "sWGS"
samples <- bind_rows(data_samples, samples)
samples[is.na(samples)] <- "HBC"
data <- data[, samples$sWGS]
data <- as.matrix(data)

### Format sample sheet
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

### Set means for each cohort
data_means <- data_freq[row.names(data), ]
data_means <- data_means[, c("all", "c2", "c3", "lof", "splice")]
colnames(data_means) <- c("Cohort", "Cluster 2", "Cluster 3", "LOF", "Splice")

data_means_mid <- as.matrix(data_means[grepl("midpoint", row.names(data_means)), ])
data_means_amp <- as.matrix(data_means[grepl("amp", row.names(data_means)), ])

### Format data into seperate features
data_cov <- as.matrix(data[grepl("coverage", row.names(data)), ])
data_mid <- as.matrix(data[grepl("midpoint", row.names(data)), ])
data_amp <- as.matrix(data[grepl("amp", row.names(data)), ])

### Set colors
col_fun <- colorRamp2(c(-6, 0, 6), 
                        c("#1f78b4", "white", "#e31a1c"))
col_means <- colorRamp2(c(-3, 0, 3), 
                        c("#1f78b4", "white", "#e31a1c"))
col_status <- c(Healthy = "#B2DF8A", Positive = "#fb9a99", Negative = "#a6cee3")
col_previous <- c(Healthy = "#B2DF8A", Yes = "#fb9a99", No = "grey95")
col_age <- c(Healthy = "#B2DF8A", Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_type <- c(Healthy = "#B2DF8A", LOF = "#FC8D62", Splice = "#66C2A5", "1" = "#E5C494", "2" = "#FFD92F", "3" = "#8DA0CB", "5" = "#E78AC3")

### Make annotations
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

right_annotation1 <- rowAnnotation("Means" = data_means_mid,
                                  show_annotation_name = TRUE,
                                  border = FALSE,
                                  col = list("Means" = col_means),
                                  annotation_name_gp = gpar(fontsize = 8),
                                  annotation_name_side = "top",
                                  annotation_name_rot = 90,
                                  simple_anno_size = unit(0.25, "cm"),
                                  show_legend = c(FALSE))

right_annotation2 <- rowAnnotation("Means" = data_means_amp,
                                   show_annotation_name = FALSE,
                                   border = FALSE,
                                   col = list("Means" = col_means),
                                   annotation_name_gp = gpar(fontsize = 8),
                                   annotation_name_side = "top",
                                   annotation_name_rot = 90,
                                   simple_anno_size = unit(0.25, "cm"),
                                   show_legend = c(FALSE))


### Make Legend
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
                                                  at = c("LOF", "Splice", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 5"),
                                                  legend_gp = gpar(fill = c("#FC8D62", "#66C2A5", "#E5C494", "#FFD92F", "#8DA0CB", "#E78AC3")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "SD from Healthy",
                                                  at = c(-6, 0, 6),
                                                  col_fun = col_fun,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm")),
                                           Legend(title = "Cohort Mean",
                                                  at = c(-3, 0, 3),
                                                  col_fun = col_means,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm"))))

### Set splits
col_split <- samples$cancer_status
col_split <- factor(col_split, levels = c("Healthy", "Negative", "Positive"))

### Plot clustering
Midpoint <- Heatmap(data_mid,
                    col = col_fun,
                    column_split = col_split,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    show_heatmap_legend = FALSE,
                    top_annotation = top_annotation,
                    right_annotation = right_annotation1,
                    column_title_gp = gpar(fontsize = 8),
                    row_labels = NULL,
                    row_title = "Midpoint Coverage",
                    row_title_gp = gpar(fontsize = 8),
                    row_title_rot = 90,
                    column_dend_height = unit(2, "cm"),
                    border = FALSE)

Amplitude <- Heatmap(data_amp,
                     col = col_fun,
                     column_split = col_split,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     show_heatmap_legend = FALSE,
                     right_annotation = right_annotation2,
                     column_title_gp = gpar(fontsize = 8),
                     row_labels = NULL,
                     row_title = "Amplitude",
                     row_title_gp = gpar(fontsize = 8),
                     row_title_rot = 90,
                     border = FALSE)

pdf(file.path(outdir, "fragment_nucleosome_TFBS_clustering_expanded.pdf"), height = 5, width = 6)
Heatmap <- Midpoint %v% Amplitude
draw(Heatmap, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "right")
dev.off()
