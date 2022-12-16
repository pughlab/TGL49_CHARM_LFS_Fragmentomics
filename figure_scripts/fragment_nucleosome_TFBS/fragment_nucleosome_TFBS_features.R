library(dplyr)
library(matrixStats)
library(ComplexHeatmap)

### Set paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/griffin/TFBS"
outdir <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_nucleosome_TFBS"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/griffin/TFBS"

data_griffin <- list.files(path, "features", full.names = TRUE)
data_samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
data_normal <- list.files(healthy_path, "features", full.names = TRUE)

### Import data 
datalist <- lapply(data_griffin, function(x){read.delim(file = x)})
data_griffin <- do.call(rbind, datalist)
datalist <- lapply(data_normal, function(x){read.delim(file = x)})
data_normal <- do.call(rbind, datalist)
data_samples <- read.delim(data_samples)

### Read in coverages and exclude <1x samples
#coverage <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/GATK_coverage/CHARM_LFS_coverage_summary.txt"
#coverage <- read.delim(coverage)
#coverage <- coverage[coverage$coverage_swgs >= 1,]
#coverage <- coverage[!(is.na(coverage$coverage_swgs)), ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
irregular <- c("TGL49_0019_Cf_U_PE_317_WG", "TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG", "TGL49_0216_Cf_n_PE_483_WG") ### As determined by PCA analysis
data_samples <- data_samples[!(data_samples$sWGS %in% exclude |
                                 data_samples$sWGS %in% irregular), ]
data_samples <- data_samples[!(data_samples$germline_mutation == "unknown"), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_griffin), ]
#data_samples <- data_samples[data_samples$sWGS %in% coverage$sWGS, ]

### Format griffin data
row.names(data_griffin) <- data_griffin$features
data_griffin <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS]
data_griffin <- as.data.frame(t(data_griffin))

row.names(data_normal) <- data_normal$features
data_normal <- data_normal[, -1]
data_normal <- as.data.frame(t(data_normal))

### combine data and run PCA
data <- bind_rows(data_griffin, data_normal)
data <- scale(data, center = TRUE, scale = TRUE)

pca <- prcomp(data, center = TRUE,scale. = TRUE)

data_pca <- as.data.frame(pca[["x"]])
pca <- data_pca
pca$diag <- ifelse(row.names(pca) %in% data_samples$sWGS, "LFS", "HBC")
pca <- merge(pca, data_samples, by.x = "row.names", by.y = "sWGS", all = TRUE)
pca[is.na(pca)] <- "HBC"

### Format PCA table
pca$cancer_status <- factor(pca$cancer_status, levels = c("HBC", "negative", "positive"),
                            labels = c("Healthy", "Negative", "Positive"))
pca$mutation_type <- factor(pca$mutation_type, levels = c("HBC", "1", "2", "3", "5", "LOF", "Splice"),
                            labels = c("Healthy", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 5", "LOF", "Splice"))

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 15), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "bottom",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 12),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Plot PCA
plot_pca <- ggplot(pca, aes(PC1, PC2, color = cancer_status)) +
  geom_point(alpha = 0.75) +
  labs(color = "Cancer Status") +
  scale_color_manual(values = c("#488f31", "grey65", "#de425b")) +
  theme
plot_pca

plot_type <- ggplot(pca, aes(PC1, PC2, color = mutation_type)) +
  geom_point(alpha = 0.75) +
  labs(color = "Mutation Type") +
  scale_color_manual(values = c("#488f31", "#FB9A99", "#A6CEE3", "#FDBF6F","#CAB2D6", "grey65", "#B2DF8A")) +
  theme
plot_type

Figure <- ggarrange(plot_pca, plot_type, align = "hv", nrow = 1)
Figure
ggsave(file.path(outdir, "fragment_nucleosome_TFBS_pca.pdf"), Figure, width = 6, height = 4)

### Cluster features
### Make a combined samples list
samples <- as.data.frame(row.names(data_normal))
colnames(samples) <- "sWGS"
samples <- bind_rows(data_samples, samples)
samples[is.na(samples)] <- "HBC"
data_pca <- t(data_pca)
data_pca <- data_pca[, samples$sWGS]

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

### Set colors
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
                                                  grid_height = unit(1, "mm"))),
                               direction = "horizontal")

### Plot clustering
pdf(file.path(outdir, "fragment_nucleosome_TFBS_clustering.pdf"), height = 2.5, width = 4)
Heatmap <- Heatmap(data_pca,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_heatmap_legend = FALSE,
                   row_order = row.names(data_pca),
                   top_annotation = top_annotation,
                   row_labels = NULL,
                   column_title_gp = gpar(fontsize = 6),
                   row_title_gp = gpar(fontsize = 0),
                   row_title_rot = 0,
                   column_dend_height = unit(2, "cm"),
                   height = unit(0, "cm"),
                   border = FALSE)
draw(Heatmap, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE, heatmap_legend_side = "bottom")
dev.off()
