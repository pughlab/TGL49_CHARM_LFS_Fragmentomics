library(gsignal)
library(dplyr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(pracma)

### Set paths
sites <- c("CTCF", "TP53", "TP53_targets_TSS", "housekeeping_TSS")
names <- c("CTCF", "TP53", "TP53 Targets TSS", "Housekeeping TSS")

path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/griffin_all"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_nucleosome_TP53"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/griffin_all"

for (i in c(1:length(sites))) {
  ### Find paths
  site <- sites[[i]]
  name <- names[[i]]
  data_griffin <- list.files(path, site, recursive = TRUE, full.names = TRUE)
  data_griffin <- data_griffin[grepl(paste0("corrected_", site, ".txt"), data_griffin)]
  data_griffin <- data_griffin[!grepl("Ulz", data_griffin)]
  data_normal <- list.files(healthy_path, site, recursive = TRUE, full.names = TRUE)
  data_normal <- data_normal[grepl(paste0("corrected_", site, ".txt"), data_normal)]
  data_normal <- data_normal[!grepl("Ulz", data_normal)]
  
  ### Import data 
  data_griffin <- read.delim(data_griffin)
  data_normal <- read.delim(data_normal)
  data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")
  
  ### Apply Savitzky-Golay filter
  distance <- data_griffin$distance
  data_griffin <- sapply(data_griffin[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_griffin <- as.data.frame(data_griffin)
  data_griffin$distance <- distance
  
  data_normal <- sapply(data_normal[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_normal <- as.data.frame(data_normal)
  data_normal$distance <- distance
  
  ### Remove failed and unknown samples and format 
  exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
  abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
  data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude, abnormal)), ]
  
  ### Adjust and center curves
  means <- colMeans(as.matrix(data_griffin))
  data_griffin <- as.data.frame(t(data_griffin))
  data_griffin <- (data_griffin - means) + 1
  data_griffin <- as.data.frame(t(data_griffin))
  
  means <- colMeans(as.matrix(data_normal))
  data_normal <- as.data.frame(t(data_normal))
  data_normal <- (data_normal - means) + 1
  data_normal <- as.data.frame(t(data_normal))
  
  ### Order based on clinical information
  data_samples <- data_samples[data_samples$sWGS %in% colnames(data_griffin), ]
  data_samples <- data_samples[order(data_samples$sample_parent,
                                     data_samples$timepoint), ]
  data_samples$diag <- "LFS"
  
  ### Format dataframes
  data_griffin <- data_griffin[, c("distance", data_samples$sWGS)]
  
  ### Calculate the healthy median
  normal_median <- rowMedians(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
  normal_sd <- rowSds(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
  
  ### Calculate the LFS medians
  c2_median <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "2"]]))
  c2_sd <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "2"]]))
  
  c3_median <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "3"]]))
  c3_sd <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "3"]]))
  
  c5_median <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "5"]]))
  c5_sd <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "5"]]))
  
  lof_median <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "LOF"]]))
  lof_sd <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "LOF"]]))
  
  sp_median <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "Splice"]]))
  sp_sd <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% data_samples$sWGS[data_samples$mutation_type == "Splice"]]))
  
  ### Make median tables
  normal_median <- as.data.frame(cbind(distance, normal_median, normal_sd))
  colnames(normal_median) <- c("distance", "median", "sd")
  normal_median$diag <- "HBC"
  
  c2_median <- as.data.frame(cbind(distance, c2_median, c2_sd))
  colnames(c2_median) <- c("distance", "median", "sd")
  c2_median$diag <- "Cluster 2"
  
  c3_median <- as.data.frame(cbind(distance, c3_median, c3_sd))
  colnames(c3_median) <- c("distance", "median", "sd")
  c3_median$diag <- "Cluster 3"
  
  c5_median <- as.data.frame(cbind(distance, c5_median, c5_sd))
  colnames(c5_median) <- c("distance", "median", "sd")
  c5_median$diag <- "Cluster 5"
  
  lof_median <- as.data.frame(cbind(distance, lof_median, lof_sd))
  colnames(lof_median) <- c("distance", "median", "sd")
  lof_median$diag <- "LOF"
  
  sp_median <- as.data.frame(cbind(distance, sp_median, sp_sd))
  colnames(sp_median) <- c("distance", "median", "sd")
  sp_median$diag <- "Splice"
  
  data_median <- bind_rows(normal_median, c2_median, c3_median, c5_median, lof_median, sp_median)
  data_median <- data_median[data_median$distance > -900 & data_median$distance < 900, ]
  data_median$diag <- factor(data_median$diag, levels = c("HBC", "LOF", "Splice", "Cluster 2", "Cluster 3", "Cluster 5"),
                             labels = c("Healthy", "LOF", "Splice", "Missense 2", "Missense 3", "Missense 5"))
  
  median <- data_median$median[data_median$diag == "Healthy" & data_median$distance == 0]

  ### Set Theme
  theme <- theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
                 axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(fill = NA),
                 panel.background = element_blank(),
                 legend.position = "none",
                 legend.key = element_rect(fill = "white"),
                 legend.text = element_text(size = 12),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 13),
                 axis.text = element_text(size = 13),
                 axis.title = element_text(size = 13),
                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ### Plot curves
  fig <- ggplot(data_median) +
    geom_line(data = data_median, aes(distance, median, group = diag), size = 0.5, color = "red") +
    geom_ribbon(data = data_median, aes(distance, ymin = median - sd, ymax = median + sd, group = diag), 
                fill = "#FB9A99", alpha = 0.5) +
    geom_hline(yintercept = median, linetype = "dashed") +
    xlab("Distance from site (bp)") + 
    ylab("Coverage") +
    ggtitle(name) + 
    facet_wrap(.~diag, nrow = 1) +
    theme +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(-900, 900), expand = c(0,0))
  fig
  
  ggsave(file.path(outdir, paste0("fragment_nucleosome_TP53_curves_", site, "_type.pdf")), fig, width = 9.5, height = 3)
}





