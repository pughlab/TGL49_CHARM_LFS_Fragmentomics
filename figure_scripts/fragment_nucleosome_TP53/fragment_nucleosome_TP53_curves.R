library(gsignal)
library(dplyr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(pracma)

### Set paths
sites <- c("CTCF", "TP53", "TP53_targets_TSS", "housekeeping_TSS", "EZH2_targets", "H3K27me3", "TP53_DMR")
names <- c("CTCF", "TP53", "TP53 Targets", "Housekeeping", "EZH2 Targets", "H3K27me3 Genes", "TP53 DMRs")

path <- "data/griffin"
outdir <- ""
healthy_path <- "hbc/griffin"

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
  data_samples <- read.delim("sample_list.txt")
  
  ### Apply Savitzky-Golay filter
  distance <- data_griffin$distance
  data_griffin <- sapply(data_griffin[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_griffin <- as.data.frame(data_griffin)
  
  data_normal <- sapply(data_normal[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_normal <- as.data.frame(data_normal)
  
  ### Remove failed and unknown samples and format 
  exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
  abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
  data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude, abnormal)), ]
  data_griffin <- data_griffin[, colnames(data_griffin) %in% c(data_samples$sWGS)]
  
  ### Adjust and center curves
  means <- colMeans(as.matrix(data_griffin))
  data_griffin <- as.data.frame(t(data_griffin))
  data_griffin <- (data_griffin - means) + 1
  data_griffin <- as.data.frame(t(data_griffin))
  data_griffin$distance <- distance
  
  means <- colMeans(as.matrix(data_normal))
  data_normal <- as.data.frame(t(data_normal))
  data_normal <- (data_normal - means) + 1
  data_normal <- as.data.frame(t(data_normal))
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
  LFS_median <- rowMedians(as.matrix(data_griffin[, !(colnames(data_griffin) == "distance")]))
  LFS_sd <- rowSds(as.matrix(data_griffin[, !(colnames(data_griffin) == "distance")]))
  
  samples_neg <- data_samples[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no", ]
  LFS_median_neg <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% samples_neg$sWGS]))
  LFS_sd_neg <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% samples_neg$sWGS]))
  
  samples_sur <- data_samples[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes", ]
  LFS_median_sur <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% samples_sur$sWGS]))
  LFS_sd_sur <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% samples_sur$sWGS]))
  
  samples_pos <- data_samples[data_samples$cancer_status == "positive", ]
  LFS_median_pos <- rowMedians(as.matrix(data_griffin[, colnames(data_griffin) %in% samples_pos$sWGS]))
  LFS_sd_pos <- rowSds(as.matrix(data_griffin[, colnames(data_griffin) %in% samples_pos$sWGS]))
  
  ### Make median tables
  data_median <- rbind(data.frame(distance = distance, median = normal_median, sd = normal_sd, diag = "Healthy"),
                       data.frame(distance = distance, median = LFS_median, sd = LFS_sd, diag = "LFS (All)"),
                       data.frame(distance = distance, median = LFS_median_neg, sd = LFS_sd_neg, diag = "LFS-H"),
                       data.frame(distance = distance, median = LFS_median_sur, sd = LFS_sd_sur, diag = "LFS-PC"),
                       data.frame(distance = distance, median = LFS_median_pos, sd = LFS_sd_pos, diag = "LFS-AC"))
  data_median$site <- name
  
  data_median <- data_median[data_median$distance > -900 & data_median$distance < 900, ]
  data_median$diag <- factor(data_median$diag, levels = c("Healthy", "LFS (All)", "LFS-H", "LFS-PC", "LFS-AC"))
  
  median <- median(data_median$median[data_median$diag == "Healthy" & data_median$distance %in% c(-30, -15, 0, 15, 30)])
  
  ### Assign variable name
  assign(paste0("data_", site), data_median)
  assign(paste0("median_", site), median)
}

### Make plotting tables
data_plot1 <- rbind(data_CTCF, data_TP53)
data_plot1$site <- factor(data_plot1$site, levels = c("CTCF", "TP53"))
data_median1 <- data.frame(median = c(median_CTCF, median_TP53),
                           site = c("CTCF", "TP53"))
data_median1$site <- factor(data_median1$site, levels = c("CTCF", "TP53"))

data_plot2 <- rbind(data_housekeeping_TSS, data_TP53_targets_TSS)
data_plot2$site <- factor(data_plot2$site, levels = c("Housekeeping", "TP53 Targets"))
data_median2 <- data.frame(median = c(median_housekeeping_TSS, median_TP53_targets_TSS),
                           site = c("Housekeeping", "TP53 Targets"))
data_median2$site <- factor(data_median2$site, levels = c("Housekeeping", "TP53 Targets"))

data_plot3 <- rbind(data_EZH2_targets, data_H3K27me3, data_TP53_DMR)
data_plot3$site <- factor(data_plot3$site, levels = c("EZH2 Targets", "H3K27me3 Genes", "TP53 DMRs"))
data_median3 <- data.frame(median = c(median_EZH2_targets, median_H3K27me3, median_TP53_DMR),
                           site = c("EZH2 Targets", "H3K27me3 Genes", "TP53 DMRs"))
data_median3$site <- factor(data_median3$site, levels = c("EZH2 Targets", "H3K27me3 Genes", "TP53 DMRs"))

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
fig1 <- ggplot(data_plot1) +
  geom_line(aes(distance, median, group = diag), size = 0.5, color = "red") +
  geom_ribbon(aes(distance, ymin = median - sd, ymax = median + sd, group = diag), 
              fill = "#FB9A99", alpha = 0.5) +
  geom_hline(data = data_median1, aes(yintercept = median), linetype = "dashed") +
  xlab("Distance from site (bp)") + 
  ylab("Coverage") +
  ggtitle("Binding Sites") + 
  facet_grid(site~diag, scales = "free_y") +
  theme +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-900, 900), expand = c(0,0))
fig1

ggsave(file.path(outdir, paste0("fragment_nucleosome_TP53_curves1.pdf")), fig1, width = 14, height = 4)

fig2 <- ggplot(data_plot2) +
  geom_line(aes(distance, median, group = diag), size = 0.5, color = "red") +
  geom_ribbon(aes(distance, ymin = median - sd, ymax = median + sd, group = diag), 
              fill = "#FB9A99", alpha = 0.5) +
  geom_hline(data = data_median2, aes(yintercept = median), linetype = "dashed") +
  xlab("Distance from site (bp)") + 
  ylab("Coverage") +
  ggtitle("Transcription Start Sites") + 
  facet_grid(site~diag, scales = "free_y") +
  theme +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-900, 900), expand = c(0,0))
fig2

ggsave(file.path(outdir, paste0("fragment_nucleosome_TP53_curves2.pdf")), fig2, width = 14, height = 4)

fig3 <- ggplot(data_plot3) +
  geom_line(aes(distance, median, group = diag), size = 0.5, color = "red") +
  geom_ribbon(aes(distance, ymin = median - sd, ymax = median + sd, group = diag), 
              fill = "#FB9A99", alpha = 0.5) +
  geom_hline(data = data_median3, aes(yintercept = median), linetype = "dashed") +
  xlab("Distance from site (bp)") + 
  ylab("Coverage") +
  ggtitle("Transcription Start Sites") + 
  facet_grid(site~diag, scales = "free_y") +
  theme +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-900, 900), expand = c(0,0))
fig3

ggsave(file.path(outdir, paste0("fragment_nucleosome_TP53_curves3.pdf")), fig3, width = 14, height = 6)



