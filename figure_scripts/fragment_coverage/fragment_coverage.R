library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(matrixStats)
library(dendextend)
library(psych)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_coverage"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"

### Import data (Starting with the 5Mb ratios)
data_cov <- read.delim(list.files(file.path(path, "fragmentomics"), "frag", full.names = TRUE))
data_normal <- read.delim(list.files(file.path(healthy, "fragmentomics"), "frag", full.names = TRUE))
samples <- read.delim(file.path(path, "samples/sample_list.txt"))
de_table <- read.delim(file.path(outdir, "fragment_coverage_differential.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]

### Keep only plasma samples and order based on clinical information
samples <- samples[!(samples$germline_mutation == "unknown"), ]
samples <- samples[samples$sWGS %in% colnames(data_cov), ]

### Set chromosomes
data_chr <- data_cov[, c("seqnames", "start", "end")]
row.names(data_chr) <- with(data_chr, paste0(seqnames, "_", start))

### Format DE table
de_table <- de_table[de_table$padj_cohort < 0.05 | 
                       de_table$padj_c2 < 0.05 | 
                       de_table$padj_c3 < 0.05 |
                       de_table$padj_lof < 0.05 |
                       de_table$padj_splice < 0.05, ]
de_table <- de_table[order(factor(de_table$location, levels = row.names(data_cov))), ]

### Format dataframe and order samples properly
row.names(data_normal) <- with(data_normal, paste0(seqnames, "_", start))
data_normal <- data_normal[, -c(1:3)]
sums <- colSums(data_normal, na.rm = TRUE)
data_normal <- sweep(data_normal, 2, sums, "/")
data_normal <- data_normal*1000000
data_normal <- data_normal[row.names(data_normal) %in% de_table$location,]
data_normal <- data_normal[complete.cases(data_normal), ]

row.names(data_cov) <- with(data_cov, paste0(seqnames, "_", start))
data_cov <- data_cov[, samples$sWGS]
sums <- colSums(data_cov, na.rm = TRUE)
data_cov <- sweep(data_cov, 2, sums, "/")
data_cov <- data_cov*1000000
data_cov <- data_cov[row.names(data_normal), ]

### Make matrix for healthy median
normal <- data.frame(median = rowMedians(as.matrix(data_normal)),
                     sd = rowSds(as.matrix(data_normal)))
row.names(normal) <- row.names(data_normal)

### Calculate the distance from the healthy median
data_change <- (data_cov - normal$median)/normal$sd
data_normal <- (data_normal - normal$median)/normal$sd

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

### Set additional annotations
top_annotation <- HeatmapAnnotation("Cancer Status" = data_status,
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
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    simple_anno_size = unit(0.25, "cm"),
                                    show_legend = c(FALSE, FALSE, FALSE, FALSE, FALSE))

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
                                                  legend_width = unit(2, "cm"))))

## Generate oncoprint
pdf(file.path(outdir, "fragment_genome_coverage.pdf"), height = 4, width = 6)
Heatmap <- Heatmap(data_change,
                   col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   top_annotation = top_annotation,
                   row_labels = NULL,
                   column_title_gp = gpar(fontsize = 6),
                   row_title_gp = gpar(fontsize = 0),
                   row_title_rot = 0,
                   row_dend_width = unit(2, "cm"),
                   cluster_row_slices = TRUE,
                   border = FALSE)
draw(Heatmap, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE)
dev.off()
