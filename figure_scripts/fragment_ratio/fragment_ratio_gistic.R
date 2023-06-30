library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(matrixStats)
library(scales)

### Set variables
path <- "data/fragment_ratio"
outdir <- ""
gistic <- "data/TCGA_merged_GISTIC_hg38.txt"

### Import data (Starting with the 5Mb ratios)
gistic <- read.delim(gistic)
de_table <- read.delim("data/fragment_ratio_differential.txt"))
data_ratio <- read.delim(file.path(path, "CHARM_LFS_ratio_5Mb.txt"))

### Set chromosomes
data_chr <- data_ratio[, c("seqnames", "arm", "start", "end")]
row.names(data_chr) <- with(data_chr, paste0(seqnames, "_", start))

### Calculate correlation to GISTIC
de_types <- de_table[, c("cohort", "c2", "c3", "lof", "splice")]
row.names(de_types) <- de_table$location

row.names(gistic) <- with(gistic, paste0(chr, "_", start))
gistic <- gistic[de_table$location, !(colnames(gistic) %in% c("chr", "arm", "start", "end"))]

types <- colnames(de_types)
cancers <- colnames(gistic)

y <- c("type", cancers)
for (type in types) {
  x <- c(type)
  for (cancer in cancers) {
    a <- corr.test(de_types[, type], gistic[, cancer])[["r"]]
    x <- c(x, a)
  }
  y <- cbind(y, x)
}

### Format correlation table
data_corr <- as.data.frame(y)
row.names(data_corr) <- data_corr$y
colnames(data_corr) <- data_corr[1, ]
data_corr <- data_corr[-1,-1]
colnames(data_corr) <- c("Cohort", "Missense 2", "Missense 3", "LOF", "Splice")
data_corr <- data_corr[colnames(gistic), ]
data_corr <- mutate_all(data_corr, function(x) as.numeric(as.character(x)))
data_corr <- as.matrix(data_corr)

### Format GISTIC data
gistic <- as.matrix(t(gistic))

### Set differential regions
de_table$regions <- ifelse(de_table$padj_cohort < 0.05, 1, 0)

### Set annotation and heatmap colours
col_heat <- colorRamp2(c(-1, 0, 1), 
                      c("#1f78b4", "white", "#e31a1c"))
col_fun <- colorRamp2(c(-2, 0, 2), 
                      c("#1f78b4", "white", "#e31a1c"))
col_fun2 <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), 
                      c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

### Set additional annotations
top_annotation <- HeatmapAnnotation("Cohort" = de_table$cohort,
                                    "Missense 2" = de_table$c2,
                                    "Missense 3" = de_table$c3,
                                    "LOF" = de_table$lof,
                                    "Splice" = de_table$splice,
                                    col = list("Cohort" = col_fun,
                                               "Missense 2" = col_fun,
                                               "Missense 3" = col_fun,
                                               "LOF" = col_fun,
                                               "Splice" = col_fun),
                                    border = FALSE,
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 8),
                                    simple_anno_size = unit(0.5, "cm"))

### Set legend labels
annotation_legend = packLegend(list = list(Legend(title = "GISTIC",
                                                  at = c(-1, 0, 1),
                                                  col_fun = col_heat,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm")),
                                           Legend(title = "SD from Healthy",
                                                  at = c(-2, 0, 2),
                                                  col_fun = col_fun,
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  direction = "horizontal",
                                                  grid_height = unit(2, "mm"),
                                                  legend_width = unit(2, "cm"))),
                               direction = "horizontal")

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
                                  "               ", "                ", "                ", "                  ","                   ","                   "))

## Set column/row order and column/row splits
chr_order <- colnames(gistic)
column_split <- data_chr$arm
row_order <- row.names(gistic)

### Make correlation heatmap
fig_corr <- Heatmap(data_corr,
                    col = col_fun2,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.2f", data_corr[i, j]), x, y, gp = gpar(fontsize = 6))
                    },
                    show_heatmap_legend = FALSE,
                    column_title_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 8),
                    row_title_gp = gpar(fontsize = 0),
                    row_names_gp = gpar(fontsize = 8),
                    height = unit(0.5*nrow(data_corr), "cm"),
                    width = unit(0.25*nrow(data_corr), "cm"))
draw(fig_corr)

## Make GISTIC Heatmap
Heatmap <- Heatmap(gistic,
                   col = col_fun2,
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   top_annotation = top_annotation,
                   column_order = chr_order,
                   #row_order = row_order,
                   column_split = column_split,
                   column_title_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8),
                   row_title_gp = gpar(fontsize = 0),
                   row_names_gp = gpar(fontsize = 8),
                   height = unit(0.5*nrow(gistic), "cm"),
                   border = FALSE)
draw(Heatmap, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "right")

### Combine heatmaps
pdf(file.path(outdir, "fragment_gistic.pdf"), height = 5, width = 10)
Figure <- Heatmap + fig_corr
draw(Figure, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "bottom")
dev.off()

