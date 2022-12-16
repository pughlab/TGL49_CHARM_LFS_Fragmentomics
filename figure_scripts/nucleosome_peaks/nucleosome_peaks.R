library(tidyverse)
library(dplyr)
library(matrixStats)
library(reshape2)
library(ggh4x)

### Set paths
sites <- c("genome", "fragments")
labels <- c("Genome", "167bp fragments")

path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/nucleosome_peaks"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/nucleosome_peaks"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/nucleosome_peaks"

for (i in c(1:length(sites))) {
  ### Set variables
  site <- sites[[i]]
  label <- labels[[i]]
  
  ### Find paths
  data <- list.files(path, ".txt", full.names = TRUE)
  data <- data[grepl(site, data)]
  data_normal <- list.files(healthy_path, ".txt", full.names = TRUE)
  data_normal <- data_normal[grepl(site, data_normal)]
  
  ### Import data 
  data <- read.delim(data)
  data_normal <- read.delim(data_normal)
  data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")
  
  ### Remove failed and unknown samples and format 
  exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
  #abnormal <- c("TGL49_0018_Cf_U_PE_323_WG")
  data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]
  
  ### Format data
  data <- data[, !(colnames(data) %in% c("chr"))]
  names <- unique(data$sample)
  data_sum <- aggregate(. ~ sample, data, sum)
  colnames(data_sum) <- c("sample", -1000:1000)
  data_sum <- data_sum[, colnames(data_sum) %in% c("sample", -300:300)]  
  data_sum[,2:ncol(data_sum)] <- data_sum[,2:ncol(data_sum)]/rowSums(data_sum[,2:ncol(data_sum)])*100
  
  data_normal <- data_normal[, !(colnames(data_normal) %in% c("chr"))]
  names <- unique(data_normal$sample)
  data_normal_sum <- aggregate(. ~ sample, data_normal, sum)
  colnames(data_normal_sum) <- c("sample", -1000:1000)
  data_normal_sum <- data_normal_sum[, colnames(data_normal_sum) %in% c("sample", -300:300)]  
  data_normal_sum[,2:ncol(data_normal_sum)] <- data_normal_sum[,2:ncol(data_normal_sum)]/rowSums(data_normal_sum[,2:ncol(data_normal_sum)])*100
  
  ### Calculate the healthy median
  normal_median <- colMedians(as.matrix(data_normal_sum[,2:ncol(data_normal_sum)]))
  normal_sd <- colSds(as.matrix(data_normal_sum[,2:ncol(data_normal_sum)]))
  
  ### Calculate the LFS median
  neg_median <- colMedians(as.matrix(data_sum[data_sum$sample %in% data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no"], 2:ncol(data_sum)]))
  neg_sd <- colSds(as.matrix(data_sum[data_sum$sample %in% data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no"], 2:ncol(data_sum)]))
  
  sur_median <- colMedians(as.matrix(data_sum[data_sum$sample %in% data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes"], 2:ncol(data_sum)]))
  sur_sd <- colSds(as.matrix(data_sum[data_sum$sample %in% data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes"], 2:ncol(data_sum)]))
  
  pos_median <- colMedians(as.matrix(data_sum[data_sum$sample %in% data_samples$sWGS[data_samples$cancer_status == "positive"], 2:ncol(data_sum)]))
  pos_sd <- colSds(as.matrix(data_sum[data_sum$sample %in% data_samples$sWGS[data_samples$cancer_status == "positive"], 2:ncol(data_sum)]))
  
  ### Calculate z-scores
  z_normal <- (normal_median - normal_median)/normal_sd
  z_neg <- (neg_median - normal_median)/normal_sd
  z_sur <- (sur_median - normal_median)/normal_sd
  z_pos <- (pos_median - normal_median)/normal_sd
  
  ### Make plotting table (medians)
  data_plot <- rbind(data.frame(length = c(-300:300),
                                score = z_neg,
                                cancer_status = "LFS-H"),
                     data.frame(length = c(-300:300),
                                score = z_sur,
                                cancer_status = "LFS-PC"),
                     data.frame(length = c(-300:300),
                                score = z_pos,
                                cancer_status = "LFS-AC"))
  
  ### Make plotting table (individuals)
  data_melt <- reshape2::melt(data_sum, id = "sample")
  data_melt$variable <- as.numeric(as.character(data_melt$variable))
  data_melt$diag <- "LFS"
  data_melt <- merge(data_melt, data_samples, by.x = "sample", by.y = "sWGS")
  data_melt$cancer_status <- ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "no", "LFS-H",
                           ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "yes", "LFS-PC", "LFS-AC"))
  data_melt$cancer_status <- factor(data_melt$cancer_status, levels =  c("LFS-H", "LFS-PC", "LFS-AC"))
  
  ### Make plotting table (difference to healthy)
  neg_dif <- neg_median - normal_median
  sur_dif <- sur_median - normal_median
  pos_dif <- pos_median - normal_median
  data_diff <- rbind(data.frame(length = c(-300:300),
                                score = neg_dif,
                                cancer_status = "LFS-H"),
                     data.frame(length = c(-300:300),
                                score = sur_dif,
                                cancer_status = "LFS-PC"),
                     data.frame(length = c(-300:300),
                                score = pos_dif,
                                cancer_status = "LFS-AC"))
  
  ### Merge plots
  data_melt <- data_melt[, c("sample", "variable", "value", "cancer_status")]
  colnames(data_melt) <- c("sample", "length", "score", "cancer_status")
  data_melt$analysis <- "Frequency"
  
  data_plot$sample <- data_plot$cancer_status
  data_plot$analysis <- "Z-score"
  
  data_diff$sample <- data_diff$cancer_status
  data_diff$analysis <- "Difference"
  
  data_plot <- bind_rows(data_plot, data_melt, data_diff)
  data_plot$analysis <- factor(data_plot$analysis, levels = c("Frequency", "Z-score", "Difference"))
  data_plot$cancer_status <- factor(data_plot$cancer_status, levels = c("LFS-H", "LFS-PC", "LFS-AC"))
  
  ### Set hlines
  data_lines <- data.frame(cancer_status = rep(c("LFS-H", "LFS-PC", "LFS-AC"), 2),
                           analysis = c(rep("Difference", 3), rep("Z-score", 3)),
                           line = c(rep(0, 6)))
  data_lines$analysis <- factor(data_lines$analysis, levels = c("Frequency", "Z-score", "Difference"))
  data_lines$cancer_status <- factor(data_lines$cancer_status, levels = c("LFS-H", "LFS-PC", "LFS-AC"))
  data_median <- data.frame(length = c(-300:300),
                            median = normal_median,
                            analysis = "Frequency")
  data_median$analysis <- factor(data_median$analysis, levels = c("Frequency", "Z-score", "Difference"))
  
  ### Set theme
  theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
                 axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(fill = NA),
                 panel.background = element_blank(),
                 legend.position = "bottom",
                 legend.key = element_rect(fill = "white"),
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 12),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 13),
                 axis.text = element_text(size = 13),
                 axis.title = element_text(size = 13))
  
  ### Plot curves
  fig <- ggplot(data_plot) +
    geom_line(aes(length, score, color = cancer_status, group = sample), size = 0.5, alpha = 0.5) +
    geom_line(data = data_median, aes(length, median), size = 0.5) +
    geom_hline(data = data_lines, aes(yintercept = line), linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(-83, 83), linetype = "dashed", size = 0.5) +
    xlab("Distance") + 
    ylab("") +
    labs(color = "") +
    facet_grid(analysis~cancer_status, scales = "free_y") +
    scale_color_manual(labels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"), 
                       values = c("Healthy" = "black", "LFS-H" = "#1F78B4", "LFS-PC" = "#33A02C", "LFS-AC" = "#E31A1C")) +
    ggtitle(paste0("Distance to Nearest Peak (", label, ")")) + 
    theme
  fig
  
  ggsave(file.path(outdir, paste0("nucleosome_peaks_", site, ".pdf")), fig, width = 8, height = 6.5)
  
  ### Repeat for healthy samples only
  ### Find paths
  data <- list.files(healthy_path, ".txt", full.names = TRUE)
  data <- data[grepl(site, data)]
  
  ### Import data 
  data <- read.delim(data)
  data_prop <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/insert_size/TGL49_HBC_proportions.txt")
  
  ### Format data
  data <- data[, !(colnames(data) %in% c("chr"))]
  names <- unique(data$sample)
  data_sum <- aggregate(. ~ sample, data, sum)
  colnames(data_sum) <- c("sample", -1000:1000)
  data_sum <- data_sum[, colnames(data_sum) %in% c("sample", -300:300)]  
  data_sum[,2:ncol(data_sum)] <- data_sum[,2:ncol(data_sum)]/rowSums(data_sum[,2:ncol(data_sum)])
  
  ### Split into low vs high healthy
  median <- median(data_prop$P.20_150.)
  samples_high <- data_prop$sample[data_prop$P.20_150. > median]
  samples_low <- data_prop$sample[data_prop$P.20_150. < median]
  
  data_high <- data_sum[data_sum$sample %in% samples_high, ]
  data_low <- data_sum[data_sum$sample %in% samples_low, ]
  
  ### Calculate the medians
  median_high <- colMedians(as.matrix(data_high[,2:ncol(data_high)]))
  sd_high <- colSds(as.matrix(data_high[,2:ncol(data_high)]))
  
  median_low <- colMedians(as.matrix(data_low[,2:ncol(data_low)]))
  sd_low <- colSds(as.matrix(data_low[,2:ncol(data_low)]))
  
  ### Make plotting table (individual curves)
  data_normal_high <- reshape2::melt(data_normal_sum[data_normal_sum$sample %in% samples_high, ], id = "sample")
  data_normal_high$diag <- "High Short Frags"
  
  data_normal_low <- reshape2::melt(data_normal_sum[data_normal_sum$sample %in% samples_low, ], id = "sample")
  data_normal_low$diag <- "Low Short Frags"
  
  data_normal_melt <- bind_rows(data_normal_high, data_normal_low)
  data_normal_melt$variable <- as.numeric(as.character(data_normal_melt$variable))
  
  ### Calculate z-scores
  z_score <- (median_high - median_low)/sd_low
  z_score <- data.frame(length = c(-300:300),
                        score = z_score,
                        analysis = "Z-score")
  
  ### Make plotting table (difference to healthy)
  data_diff <- median_high - median_low
  data_diff <- data.frame(length = c(-300:300),
                          score = data_diff,
                          analysis = "Difference")
  
  ### Merge plots
  colnames(data_normal_melt) <- c("sample", "length", "score", "cancer_status")
  data_normal_melt$analysis <- "Frequency"
  
  z_score$sample <- "Low Short Fragments"
  z_score$cancer_status <- z_score$sample
  
  data_diff$sample <- "Low Short Fragments"
  data_diff$cancer_status <- data_diff$sample
  
  data_plot <- bind_rows(z_score, data_normal_melt, data_diff)
  data_plot$analysis <- factor(data_plot$analysis, levels = c("Frequency", "Z-score", "Difference"))
  data_plot$cancer_status <- factor(data_plot$cancer_status, levels = c("Low Short Frags", "High Short Frags"))
  
  ### Plot curves
  fig <- ggplot(data_plot) +
    geom_line(aes(length, score, color = cancer_status), size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(-80, 83), linetype = "dashed", size = 0.5) +
    xlab("Distance") + 
    ylab("") +
    labs(color = "") +
    facet_grid(analysis~., scales = "free_y") +
    facetted_pos_scales(y = list(NULL, NULL, scale_y_continuous(limits = c(-0.01, 0.005)))) +
    scale_color_manual(values = c("black", "#E31A1C")) +
    ggtitle(paste0("Healthy Controls (", label, ")")) + 
    theme
  fig
  
  ggsave(file.path(outdir, paste0("nucleosome_peaks_hbc_", site, ".pdf")), fig, width = 4, height = 6.5)
}


