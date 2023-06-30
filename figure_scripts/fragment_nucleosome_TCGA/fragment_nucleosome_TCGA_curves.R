library(gsignal)
library(dplyr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(pracma)

### Set paths
sites <- c("PRAD", "BRCA", "BLCA", "LUAD")
names <- c("Prostate Cancer", "Breast Cancer", "Bladder Cancer", "Lung Cancer")
lists <- c("prostate", "breast", "bladder", "lung")

path <- "data/griffin/TCGA"
outdir <- ""
healthy_path <- "hbc/griffin/TCGA"
samples_path <- "sample_list.txt"

for (i in c(1:length(sites))) {
  ### Set variables
  site <- sites[[i]]
  name <- names[[i]]
  list <- lists[[i]]
  
  ### Find paths
  data_site <- list.files(path, site, full.names = TRUE)
  data_site <- data_site[grepl("corrected", data_site)]
  normal_site <- list.files(healthy_path, site, full.names = TRUE)
  normal_site <- normal_site[grepl("corrected", normal_site)]
  
  ### Import data 
  data_griffin <- read.delim(data_site)
  data_normal <- read.delim(normal_site)
  data_samples <- read.delim(samples_path)
  
  ### Apply Savitzky-Golay filter
  distance <- data_griffin$distance
  data_griffin <- sapply(data_griffin[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_griffin <- as.data.frame(data_griffin)
  
  data_normal <- sapply(data_normal[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_normal <- as.data.frame(data_normal)
  
  ### Remove failed and unknown samples and format 
  exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
  abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG", "TGL49_0019_Cf_U_PE_317_WG")
  data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude, abnormal)), ]
  
  ### Adjust and center curves
  #means <- colMeans(as.matrix(data_griffin))
  #data_griffin <- as.data.frame(t(data_griffin))
  #data_griffin <- (data_griffin - means) + 1
  #data_griffin <- as.data.frame(t(data_griffin))
  means <- colMeans(data_griffin[c(54,55,56,77,78,79), ])
  means <- means - mean(means)
  data_griffin <- sweep(data_griffin, 2, means)
  data_griffin$distance <- distance
  
  means <- colMeans(data_normal[c(54,55,56,77,78,79), ])
  means <- means - mean(means)
  data_normal <- sweep(data_normal, 2, means)
  data_normal$distance <- distance
  
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
  
  ### Calculate the LFS median
  neg_samples <- data_samples[data_samples$cancer_status == "negative", ]
  neg_griffin <- data_griffin[, c("distance", neg_samples$sWGS)]
  LFS_median <- rowMedians(as.matrix(neg_griffin[, !(colnames(neg_griffin) == "distance")]))
  LFS_sd <- rowSds(as.matrix(neg_griffin[, !(colnames(neg_griffin) == "distance")]))
  
  ### Make median tables
  normal_median <- as.data.frame(cbind(distance, normal_median, normal_sd))
  colnames(normal_median) <- c("distance", "median", "sd")
  normal_median$diag <- "Healthy"
  
  LFS_median <- as.data.frame(cbind(distance, LFS_median, LFS_sd))
  colnames(LFS_median) <- c("distance", "median", "sd")
  LFS_median$diag <- "LFS Negative"
  
  data_median <- bind_rows(LFS_median, normal_median)
  
  ### Merge and melt dataframes
  samples <- data_samples[data_samples$cancer_type == list, ]
  data <- data_griffin[, c("distance", samples$sWGS)]
  #data$distance <- data$distance - 8.5
  data_melt <- reshape2::melt(data, id = "distance")
  data_melt <- merge(data_melt, samples, by.x = "variable", by.y = "sWGS", all = TRUE)
  data_melt <- data_melt[data_melt$distance > -901 & data_melt$distance < 901, ]
  data_melt$diag <- "LFS Negative"
  
  ### Set median line
  median <- min(data_median$median[data_median$diag == "Healthy" & data_median$distance %in% c("-30", "-15", "0", "15", "30")])
  
  ### Set Theme
  theme <- theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), 
                 axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 legend.position = "right",
                 legend.key = element_rect(fill = "white"),
                 legend.text = element_text(size = 12),
                 legend.title = element_text(size = 12),
                 strip.background = element_rect(fill = "white"),
                 strip.text = element_text(size = 13),
                 axis.text = element_text(size = 13),
                 axis.title = element_text(size = 13),
                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ### Plot curves
  fig <- ggplot(data_median) +
    geom_line(aes(distance, median), color = "black", size = 0.5) +
    geom_ribbon(aes(distance, ymin = median - sd, ymax = median + sd), fill = "black", alpha = 0.5) +
    geom_line(data = data_melt, aes(distance, value, color = stage, group = variable), alpha = 0.5, size = 0.5) +
    geom_hline(yintercept = median, linetype = "dashed", size = 0.5) +
    ggtitle(name) +
    xlab("Distance from site (bp)") + 
    ylab("Coverage") +
    labs(color = "LFS Cancer Positive", fill = "") +
    guides(alpha = "none") +
    scale_color_manual(labels = c("low" = "Stage 0/I/II", "high" = "Stage III/IV"),
                       values = c("low" = "#FF9333", "high" = "#339FFF")) +
    scale_alpha_manual(values = c("negative" = 0.25, "postiive" = 0.75)) +
    facet_wrap(.~diag) +
    theme +
    scale_x_continuous(limits = c(-900, 900), expand = c(0,0))
  fig
  
  ggsave(file.path(outdir, paste0("fragment_nucleosome_curves_TCGA_", site, ".pdf")), fig, width = 6, height = 4)
  plot_name <- paste0(site, "_plot")
  assign(plot_name, fig)
}

Figure <- ggarrange(PRAD_plot, BLCA_plot, nrow = 2, common.legend = TRUE, legend = "bottom")
Figure
ggsave(file.path(outdir, paste0("fragment_nucleosome_curves_TCGA_figure.pdf")), Figure, width = 5, height = 6)
