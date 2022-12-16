library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(matrixStats)
library(scales)
library(psych)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_ratio"
gistic <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/External_data/TCGA_GISTIC/TCGA_merged_GISTIC_hg38.txt"

### Import data (Starting with the 5Mb ratios)
gistic <- read.delim(gistic)
data_ratio <- read.delim(file.path(path, "fragmentomics", "CHARM_LFS_ratio_5Mb.txt"))
data_ichor <- read.delim(file.path(path, "ichorCNA", "CHARM_LFS_ichorCNA_summary_reviewed.txt"))
de_table <- read.delim(file.path(outdir, "fragment_ratio_differential.txt"))
samples <- read.delim(file.path(path, "samples/sample_list.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]

### Keep only plasma samples and order based on clinical information
samples <- samples[samples$sWGS %in% colnames(data_ratio), ]

### Set chromosomes
data_chr <- data_ratio[, c("seqnames", "arm", "start", "end")]
row.names(data_chr) <- with(data_chr, paste0(seqnames, "_", start))

### Format ratio sheet and order samples properly
row.names(data_ratio) <- with(data_ratio, paste0(seqnames, "_", start))
data_ratio <- data_ratio[, samples$sWGS]

### Make matrix for healthy median
normal <- de_table[, c("healthy_median", "healthy_sd")]
colnames(normal) <- c("median", "sd")
row.names(normal) <- de_table$location

### Calculate the distance from the healthy median
data_change <- (data_ratio - normal$median)/normal$sd

## Keep high ctDNA fraction cancer positive samples
samples <- samples[samples$cancer_status == "positive", ]
samples <- samples[samples$cancer_type %in% c("adrenal", "breast", "prostate", "bladder", "lung"), ]
samples <- merge(samples, data_ichor[, 1:5], by.x = "sWGS", by.y = "sample")
samples$TF_sum <- ifelse(samples$TF > samples$TF_short, samples$TF, samples$TF_short)
samples$cancer_type <- factor(samples$cancer_type, levels = c("adrenal", "breast", "bladder", "lung", "prostate"),
                              labels = c("ACC", "BRCA", "BLCA", "LUAD", "PRAD"))
samples <- samples[order(samples$cancer_type, 
                         samples$TF_sum), ]

data_change <- data_change[, samples$sWGS]

### Calculate correlation to GISTIC
row.names(gistic) <- with(gistic, paste0(chr, "_", start))
gistic <- gistic[de_table$location, !(colnames(gistic) %in% c("chr", "arm", "start", "end"))]

names <- colnames(data_change)
cancers <- colnames(gistic)

y <- c("type", cancers)
for (name in names) {
  x <- c(name)
  for (cancer in cancers) {
    a <- corr.test(data_change[, name], gistic[, cancer])[["r"]]
    x <- c(x, a)
  }
  y <- cbind(y, x)
}

### Format correlation table
data_corr <- as.data.frame(y)
row.names(data_corr) <- data_corr$y
colnames(data_corr) <- data_corr[1, ]
data_corr <- data_corr[-1,-1]
data_corr <- mutate_all(data_corr, function(x) as.numeric(as.character(x)))
data_corr <- as.matrix(t(data_corr))

### Set annotation and heatmap colours
col_heat <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), 
                       c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
col_type <- c(ACC = "#7eb0d5", BLCA = "#b2e061", BRCA = "#bd7ebe", LUAD = "#fdcce5", PRAD = "#b3d4ff")

### Set additional annotations
left_annotation <- rowAnnotation("Cancer Type" = samples$cancer_type,
                                 col = list("Cancer Type" = col_type),
                                 border = FALSE,
                                 show_annotation_name = TRUE,
                                 annotation_name_rot = 0,
                                 annotation_name_gp = gpar(fontsize = 0),
                                 simple_anno_size = unit(0.5, "cm"),
                                 show_legend = c(FALSE))

right_annotation <- rowAnnotation("Tumour Fraction" = anno_points(samples$TF_sum,
                                                                  pch = 16,
                                                                  size = unit(1, "mm"),
                                                                  ylim = c(0, 0.5),
                                                                  border = TRUE,
                                                                  axis_param = list(side = "bottom",
                                                                                    labels_rot = 0,
                                                                                    gp = gpar(fontsize = 6))),
                                  width = unit(1.5, "cm"),
                                  show_annotation_name = TRUE,
                                  border = TRUE,
                                  annotation_name_gp = gpar(fontsize = 8),
                                  annotation_name_side = "top",
                                  annotation_name_rot = 0)

### Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "Correlation",
                                                  at = c(-1, 0, 1),
                                                  col_fun = col_heat,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm")),
                                           Legend(title = "Cancer Type", 
                                                  at = c("ACC", "BLCA", "BRCA", "LUAD", "PRAD"),
                                                  legend_gp = gpar(fill = c("#7eb0d5", "#b2e061", "#bd7ebe", "#fdcce5", "#b3d4ff")),
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  grid_height = unit(1, "mm"),
                                                  nrow = 2)),
                               direction = "horizontal")


## Set column/row order and column/row splits
col_order <- colnames(data_corr)
row_order <- row.names(data_corr)

### Make correlation heatmap
pdf(file.path(outdir, "fragment_gistic_positive.pdf"), height = 5, width = 6)
fig_corr <- Heatmap(data_corr,
                    col = col_heat,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.2f", data_corr[i, j]), x, y, gp = gpar(fontsize = 6))
                    },
                    row_order = row_order,
                    column_order = col_order,
                    show_heatmap_legend = FALSE,
                    left_annotation = left_annotation,
                    right_annotation = right_annotation,
                    show_row_names = FALSE,
                    column_title_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 8),
                    row_title_gp = gpar(fontsize = 0),
                    row_names_gp = gpar(fontsize = 8),
                    height = unit(0.4*nrow(data_corr), "cm"),
                    width = unit(0.5*nrow(data_corr), "cm"))

draw(fig_corr, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "bottom")
dev.off()

