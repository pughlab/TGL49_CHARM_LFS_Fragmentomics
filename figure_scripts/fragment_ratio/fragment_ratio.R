library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(matrixStats)
library(dendextend)
library(psych)

### Set variables
path <- "data/fragment_ratio"
outdir <- ""
healthy <- "hbc/fragment_ratio"

### Import data (Starting with the 5Mb ratios)
data_ratio <- read.delim(file.path(path, "CHARM_LFS_ratio_5Mb.txt"))
data_normal <- read.delim(file.path(healthy, "TGL49_HBC_ratio_5Mb.txt"))
samples <- read.delim("sample_list.txt"))
de_table <- read.delim("data/fragment_ratio_differential.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]

### Keep only plasma samples and order based on clinical information
samples <- samples[!(samples$germline_mutation == "unknown"), ]
samples <- samples[samples$sWGS %in% colnames(data_ratio), ]

### Set chromosomes
data_chr <- data_ratio[, c("seqnames", "arm", "start", "end")]
row.names(data_chr) <- with(data_chr, paste0(seqnames, "_", start))

### Format ratio sheet and order samples properly
row.names(data_ratio) <- with(data_ratio, paste0(seqnames, "_", start))
data_ratio <- data_ratio[, samples$sWGS]

row.names(data_normal) <- with(data_normal, paste0(seqnames, "_", start))
data_normal <- as.matrix(data_normal[, -c(1:4)])

### Format DE table
de_table <- de_table[order(factor(de_table$location, levels = row.names(data_ratio))), ]

### Make matrix for healthy median
normal <- de_table[, c("healthy_median", "healthy_sd")]
colnames(normal) <- c("median", "sd")
row.names(normal) <- de_table$location

### Calculate the correlation between healthy median and samples
corr_normal <- corr.test(data_normal, normal$median)[["r"]]
data_corr <- corr.test(as.matrix(data_ratio), normal$median)[["r"]]
data_corr <- rbind(data_corr, corr_normal)

### Calculate the distance from the healthy median
data_change <- (data_ratio - normal$median)/normal$sd
data_normal <- (data_normal - normal$median)/normal$sd

### Set Median change of LFS cohort from HBCs
data_freq <- de_table$cohort
data_freq <- as.matrix(data_freq)
row.names(data_freq) <- row.names(data_change)

data_freq2 <- data_change
data_freq2 <- ifelse(data_freq2 > 3 | data_freq2 < -3, 1, 0)
data_freq2 <- (rowSums2(as.matrix(data_freq2), na.rm = TRUE)/ncol(data_freq2))*100
data_freq2 <- as.matrix(data_freq2)
row.names(data_freq2) <- row.names(data_change)

data_freq3 <- as.matrix(de_table$padj_cohort)
row.names(data_freq3) <- de_table$location
data_freq3 <- ifelse(data_freq3 < 0.05, 1, 0)

### Combine dataframes
healthy_samples <- as.data.frame(colnames(data_normal))
colnames(healthy_samples) <- c("sWGS")
healthy_samples$Age <- "HBC"
healthy_samples$cancer_status <- "HBC"
healthy_samples$previous_cancer <- "HBC"
healthy_samples$mutation_type <- "HBC"
healthy_samples$stage <- "HBC"

samples <- bind_rows(samples, healthy_samples)

order <- row.names(data_change)
data_change <- merge(data_change, data_normal, by = "row.names")
row.names(data_change) <- data_change$Row.names
data_change <- data_change[order, samples$sWGS]

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
samples$stage <- factor(samples$stage, levels = c("HBC", "high", "low", "unknown", ""),
                        labels = c("Healthy", "Stage III/IV", "Stage 0/I/II", "Unknown", "none"))

### Set clinical information
data_status <- as.matrix(samples$cancer_status)
row.names(data_status) <- samples$sWGS

data_age <- as.matrix(samples$Age)
row.names(data_age) <- samples$sWGS

data_type <- as.matrix(samples$mutation_type)
row.names(data_type) <- samples$sWGS

data_previous <- as.matrix(samples$previous_cancer)
row.names(data_previous) <- samples$sWGS

data_stage <- as.matrix(samples$stage)
row.names(data_stage) <- samples$sWGS

normal_median <- as.matrix(normal$median)
row.names(normal_median) <- row.names(data_change)

### Set annotation and heatmap colours
col_heat <- colorRamp2(c(-6, -3, 0, 3, 6), 
                      c("#1f78b4", "white", "white", "white", "#e31a1c"))
col_fun <- colorRamp2(c(-3, 0, 3), 
                      c("#1f78b4", "white", "#e31a1c"))
col_fun2 <- colorRamp2(c(0, 50), 
                       c("white", "#e31a1c"))
col_fun3 <- colorRamp2(c(0, 1), 
                       c("white", "#e31a1c"))
col_status <- c(Healthy = "#B2DF8A", Positive = "#fb9a99", Negative = "#a6cee3")
col_previous <- c(Healthy = "#B2DF8A", Yes = "#fb9a99", No = "grey95")
col_age <- c(Healthy = "#B2DF8A", Adult = "#6A3D9A", Pediatric = "#CAB2D6")
col_stage <- c(Healthy = "#B2DF8A", "Stage III/IV" = "#fb9a99", "Stage 0/I/II" = "#ebdc78", "none" = "grey95", Unknown = "grey65")
col_type <- c(Healthy = "#B2DF8A", LOF = "#FC8D62", Splice = "#66C2A5", "1" = "#E5C494", "2" = "#FFD92F", "3" = "#8DA0CB", "5" = "#E78AC3")
col_cor <- ifelse(samples$cancer_status == "Positive" & samples$stage == "Stage III/IV", "red", "black")

### Set additional annotations
left_annotation <- rowAnnotation("Cancer Status" = data_status,
                                 "Patient Type" = data_age,
                                 "Cancer History" = data_previous,
                                 "Tumor Stage" = data_stage,
                                 "Mutation Type" = data_type,
                                 show_annotation_name = TRUE,
                                 border = FALSE,
                                 col = list("Cancer Status" = col_status, 
                                            "Patient Type" = col_age,
                                            "Cancer History" = col_previous,
                                            "Mutation Type" = col_type,
                                            "Tumor Stage" = col_stage),
                                 annotation_name_gp = gpar(fontsize = 8),
                                 annotation_name_side = "top",
                                 annotation_name_rot = 90,
                                 simple_anno_size = unit(0.25, "cm"),
                                 show_legend = c(FALSE, FALSE, FALSE, FALSE, FALSE))

right_annotation <- rowAnnotation("Correlation to Healthy" = anno_points(data_corr,
                                                                         pch = 16,
                                                                         gp = gpar(col = col_cor),
                                                                         size = unit(0.75, "mm"),
                                                                         ylim = c(0, 1),
                                                                         border = FALSE,
                                                                         axis_param = list(side = "bottom",
                                                                                           labels_rot = 0,
                                                                                           gp = gpar(fontsize = 6))),
                                    width = unit(2, "cm"),
                                    show_annotation_name = TRUE,
                                    border = FALSE,
                                    annotation_name_gp = gpar(fontsize = 8),
                                    annotation_name_side = "bottom",
                                    annotation_name_rot = 0)

top_annotation <- HeatmapAnnotation("Healthy Median" = anno_lines(normal_median,
                                                          ylim = c(-0.08, 0.08),
                                                          axis_param = list(at = c(-0.05, 0, 0.05),
                                                                            side = "right",
                                                                            labels_rot = 0,
                                                                            gp = gpar(fontsize = 6)),
                                                          border = FALSE),
                                    "Differential Regions" = data_freq3,
                                    "% LFS >3 SD" = data_freq2,
                                    "LFS Cohort Median" = data_freq,
                                    col = list("LFS Cohort Median" = col_fun,
                                               "% LFS >3 SD" = col_fun2,
                                               "Differential Regions" = col_fun3),
                                    border = FALSE,
                                    height = unit(1.75, "cm"),
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "right",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 8),
                                    simple_anno_size = unit(0.25, "cm"))

### Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Cancer Status",
                                                  at = c("Healthy", "Negative", "Positive"),
                                                  legend_gp = gpar(fill = c("#B2DF8A", "#a6cee3", "#fb9a99")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Patient Type", 
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
                                           Legend(title = "Tumor Stage",
                                                  at = c("Stage 0/I/II", "Stage III/IV", "Unknown/None"),
                                                  legend_gp = gpar(fill = c("#ebdc78", "#fb9a99", "grey65")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "Mutation Type",
                                                  at = c("LOF", "Splice", "Missense 1", "Missense 2", "Missense 3", "Missense 5"),
                                                  legend_gp = gpar(fill = c("#FC8D62", "#66C2A5", "#E5C494", "#FFD92F", "#8DA0CB", "#E78AC3")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm")),
                                           Legend(title = "SD from Healthy\n(Heatmap)",
                                                  at = c(-6, -3, 0, 3, 6),
                                                  col_fun = col_heat,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm")),
                                           Legend(title = "SD from Healthy\n(LFS Median)",
                                                  at = c(-3, 0, 3),
                                                  col_fun = col_fun,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm")),
                                           Legend(title = "% of Cohort",
                                                  at = c(0, 50),
                                                  col_fun = col_fun2,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm"))))

## Set chromosome order
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
data_chr$arm <- factor(data_chr$arm, levels=armlevels)
data_chr$arm <- factor(data_chr$arm, levels=armlevels,
                       labels = c("","1"," ","2","  ","3","   ","4","    ","5","     ","6",
                                  "      ","7","       ","8", "        ", "9","         ","10","          ","11", "           ",
                                  "12","13","14","15","            ","16", "             ","17","              ","18",
                                  "               ", "19", "                ", "20","                   ","22"))

## Set column/row order and column/row splits
chr_order <- row.names(data_change)
column_split <- data_chr$arm

row_split <- samples$cancer_status
row_split <- factor(row_split, levels = c("Positive", "Negative", "Healthy"))

### Transpose data for heatmap
data_change <- t(data_change)

## Generate oncoprint
pdf(file.path(outdir, "fragment_genome_1.pdf"), height = 6, width = 12)
Heatmap <- Heatmap(data_change,
                   col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   top_annotation = top_annotation,
                   left_annotation = left_annotation,
                   right_annotation = right_annotation,
                   column_order = chr_order,
                   row_labels = NULL,
                   column_split = column_split,
                   row_split = row_split,
                   column_title_gp = gpar(fontsize = 6),
                   row_title_gp = gpar(fontsize = 0),
                   row_title_rot = 0,
                   row_dend_width = unit(2, "cm"),
                   cluster_row_slices = TRUE,
                   border = FALSE)
draw(Heatmap, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE)
dev.off()

### Print out correlation scores
data_corr <- as.data.frame(data_corr)
data_corr$sample <- row.names(data_corr)
colnames(data_corr) <- c("correlation", "sample")
write.table(data_corr, file.path(outdir, "fragment_genome_correlation.txt"), sep = "\t", row.names = FALSE)

