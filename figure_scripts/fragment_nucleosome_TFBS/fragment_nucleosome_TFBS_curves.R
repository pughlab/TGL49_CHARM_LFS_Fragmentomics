library(gsignal)
library(dplyr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(pracma)

### Set paths
sites <- c("TP53", "GRHL2", "ETV2", "POU5F1")

path <- "data/griffin"
outdir <- ""
healthy_path <- "hbc/griffin_all"

for (i in c(1:length(sites))) {
  ### Find paths
  site <- sites[[i]]
  data_griffin <- list.files(path, site, recursive = TRUE, full.names = TRUE)
  data_griffin <- data_griffin[grepl(paste0("corrected_", site, ".txt"), data_griffin)]
  data_normal <- list.files(healthy_path, site, recursive = TRUE, full.names = TRUE)
  data_normal <- data_normal[grepl(paste0("corrected_", site, ".txt"), data_normal)]
  
  ### Import data 
  data_griffin <- read.delim(data_griffin)
  data_normal <- read.delim(data_normal)
  data_samples <- read.delim("sample_list.txt")
  
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
  
  ### Format dataframes
  data_griffin <- data_griffin[, c("distance", data_samples$sWGS)]
  
  ### Calculate the healthy median
  normal_median <- rowMedians(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
  normal_sd <- rowSds(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
  
  ### Calculate the LFS medians
  samples_neg <- data_samples[data_samples$cancer_status == "negative", ]
  neg_median <- rowMedians(as.matrix(data_griffin[, samples_neg$sWGS]))
  neg_sd <- rowSds(as.matrix(data_griffin[, samples_neg$sWGS]))
  
  samples_pos <- data_samples[data_samples$cancer_status == "positive", ]
  pos_median <- rowMedians(as.matrix(data_griffin[, samples_pos$sWGS]))
  pos_sd <- rowSds(as.matrix(data_griffin[, samples_pos$sWGS]))
  
  ### Make median tables
  normal_median <- as.data.frame(cbind(distance, normal_median, normal_sd))
  colnames(normal_median) <- c("distance", "median", "sd")
  normal_median$diag <- "Healthy"
  
  neg_median <- as.data.frame(cbind(distance, neg_median, neg_sd))
  colnames(neg_median) <- c("distance", "median", "sd")
  neg_median$diag <- "LFS Negative"
  
  pos_median <- as.data.frame(cbind(distance, pos_median, pos_sd))
  colnames(pos_median) <- c("distance", "median", "sd")
  pos_median$diag <- "LFS Positive"
  
  data_median <- bind_rows(neg_median, pos_median, normal_median)
  data_median <- data_median[data_median$distance > -900 & data_median$distance < 900, ]
  
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
    ggtitle(site) + 
    facet_wrap(.~diag) +
    theme +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(limits = c(-900, 900), expand = c(0,0))
  fig
  
  ggsave(file.path(outdir, paste0("fragment_nucleosome_curves_", site, ".pdf")), fig, width = 9, height = 3)
}





