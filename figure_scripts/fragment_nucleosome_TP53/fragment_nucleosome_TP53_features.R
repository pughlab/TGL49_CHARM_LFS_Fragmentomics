library(dplyr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(ggpubr)
library(pracma)

### Set paths
sites <- c("CTCF", "TP53", "TP53_targets_TSS", "housekeeping_TSS")
names <- c("CTCF", "TP53", "TP53 Targets TSS", "Housekeeping TSS")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/griffin_all"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
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
  data_samples <- read.delim(samples)
  
  ### Apply Savitzky-Golay filter
  distance <- data_griffin$distance
  data_griffin <- sapply(data_griffin[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_griffin <- as.data.frame(data_griffin)
  data_griffin$distance <- distance
  
  data_normal <- sapply(data_normal[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_normal <- as.data.frame(data_normal)
  data_normal$distance <- distance
  
  ### Adjust and center curves
  means <- colMeans(as.matrix(data_griffin))
  data_griffin <- as.data.frame(t(data_griffin))
  data_griffin <- (data_griffin - means) + 1
  data_griffin <- as.data.frame(t(data_griffin))
  
  means <- colMeans(as.matrix(data_normal))
  data_normal <- as.data.frame(t(data_normal))
  data_normal <- (data_normal - means) + 1
  data_normal <- as.data.frame(t(data_normal))
  
  ### Remove failed and unknown samples and format 
  exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
  abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG", "TGL49_0019_Cf_U_PE_317_WG")
  data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude, abnormal)), ]
  data_samples <- data_samples[!(data_samples$germline_mutation == "unknown"), ]
  
  ### Order based on clinical information
  data_samples <- data_samples[data_samples$sWGS %in% colnames(data_griffin), ]
  data_samples <- data_samples[order(data_samples$sample_parent,
                                     data_samples$timepoint), ]
  
  ### Format dataframes
  row.names(data_griffin) <- data_griffin$distance
  data_griffin <- data_griffin[, data_samples$sWGS]
  data_griffin <- as.matrix(data_griffin)
  
  row.names(data_normal) <- data_normal$distance
  data_normal <- as.matrix(data_normal[, !(colnames(data_normal) == "distance")])
  
  row.names(data_griffin) <- distance
  row.names(data_normal) <- distance
  
  ### Calculate the mean central coverage
  normal_mid <- colMeans2(data_normal[row.names(data_normal) %in% c("-30", "-15", "0", "15", "30"), ])
  lfs_mid <- colMeans2(data_griffin[row.names(data_griffin) %in% c("-30", "-15", "0", "15", "30"), ])
  
  ### Calculate the amplitude
  normal_amp <- GeneCycle::periodogram(data_normal)[["spec"]]
  normal_amp <- colMaxs(normal_amp)
  lfs_amp <- GeneCycle::periodogram(data_griffin)[["spec"]]
  lfs_amp <- colMaxs(lfs_amp)
  
  ### Make table of scores
  data_scores <- as.data.frame(cbind(lfs_mid, lfs_amp))
  colnames(data_scores) <- c("mid", "amp")
  data_scores$sample <- colnames(data_griffin)
  
  data_scores_normal <- as.data.frame(cbind(normal_mid, normal_amp))
  colnames(data_scores_normal) <- c("mid", "amp")
  data_scores_normal$sample <- colnames(data_normal)
  
  data_scores <- bind_rows(data_scores, data_scores_normal)
  
  ### Remove outliers
  data_scores$mid <- ifelse(data_scores$mid > mean(data_scores$mid) + 4*sd(data_scores$mid) |
                              data_scores$mid < mean(data_scores$mid) - 4*sd(data_scores$mid), NA, data_scores$mid)
  data_scores$amp <- ifelse(data_scores$amp > mean(data_scores$amp) + 4*sd(data_scores$amp) |
                              data_scores$amp < mean(data_scores$amp) - 4*sd(data_scores$amp), NA, data_scores$amp)
  
  ### Combine sample dataframes and merge with scores
  healthy_samples <- as.data.frame(colnames(data_normal))
  colnames(healthy_samples) <- "sWGS"
  healthy_samples$diag <- "HBC"
  
  data_samples$diag <- "LFS"
  data_samples <- bind_rows(data_samples, healthy_samples)
  data_samples[is.na(data_samples)] <- "HBC"
  data_samples <- data_samples[data_samples$sWGS %in% data_scores$sample, ]
  
  data_scores <- merge(data_scores, data_samples, by.x = "sample", by.y = "sWGS")
  data_scores$mutation_type <- factor(data_scores$mutation_type, levels = c("HBC", "1", "2", "3", "5", "LOF", "Splice"),
                                      labels = c("Healthy", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 5", "LOF", "Splice"))
  
  ### Calculate statistics
  # midpoint coverage
  data_stats_mid <- data_scores %>%
    group_by(diag)%>% 
    dplyr::summarise(Median=median(mid, na.rm = TRUE),
                     mean=mean(mid, na.rm = TRUE),
                     SD=sd(mid, na.rm = TRUE),
                     N=n())
  a <- t.test(data_scores$mid[data_scores$diag=="HBC"], data_scores$mid[data_scores$diag=="LFS"])$p.value
  t_test <- c(1, a)
  data_stats_mid$pvalue <- t_test
  data_stats_mid$annot <- ifelse(data_stats_mid$pvalue < 0.05 & data_stats_mid$pvalue > 0.01, "*",
                                 ifelse(data_stats_mid$pvalue < 0.01 & data_stats_mid$pvalue > 0.001, "**",
                                        ifelse(data_stats_mid$pvalue < 0.001, "***", "")))
  
  # amplitude
  data_stats_amp <- data_scores %>%
    group_by(diag)%>% 
    dplyr::summarise(Median=median(amp, na.rm = TRUE),
                     mean=mean(amp, na.rm = TRUE),
                     SD=sd(amp, na.rm = TRUE),
                     N=n())
  a <- t.test(data_scores$amp[data_scores$diag=="HBC"], data_scores$amp[data_scores$diag=="LFS"])$p.value
  t_test <- c(1, a)
  data_stats_amp$pvalue <- t_test
  data_stats_amp$annot <- ifelse(data_stats_amp$pvalue < 0.05 & data_stats_amp$pvalue > 0.01, "*",
                                 ifelse(data_stats_amp$pvalue < 0.01 & data_stats_amp$pvalue > 0.001, "**",
                                        ifelse(data_stats_amp$pvalue < 0.001, "***", "")))
  
  ### Calculate statistics (type)
  types <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 5", "LOF", "Splice")
  
  # midpoint coverage
  data_stats_mid_type <- data_scores %>%
    group_by(mutation_type)%>% 
    dplyr::summarise(Median=median(mid, na.rm = TRUE),
                     mean=mean(mid, na.rm = TRUE),
                     SD=sd(mid, na.rm = TRUE),
                     N=n())
  t_test <- c(1)
  for (type in types) {
    a <- t.test(data_scores$mid[data_scores$mutation_type=="Healthy"], data_scores$mid[data_scores$mutation_type==type])$p.value
    t_test <- c(t_test, a)
  }
  t_test <- p.adjust(t_test)
  data_stats_mid_type$pvalue <- t_test
  data_stats_mid_type$annot <- ifelse(data_stats_mid_type$pvalue < 0.05 & data_stats_mid_type$pvalue > 0.01, "*",
                                 ifelse(data_stats_mid_type$pvalue < 0.01 & data_stats_mid_type$pvalue > 0.001, "**",
                                        ifelse(data_stats_mid_type$pvalue < 0.001, "***", "")))
  
  # amplitude
  data_stats_amp_type <- data_scores %>%
    group_by(mutation_type)%>% 
    dplyr::summarise(Median=median(amp, na.rm = TRUE),
                     mean=mean(amp, na.rm = TRUE),
                     SD=sd(amp, na.rm = TRUE),
                     N=n())
  t_test <- c(1)
  for (type in types) {
    a <- t.test(data_scores$amp[data_scores$mutation_type=="Healthy"], data_scores$amp[data_scores$mutation_type==type])$p.value
    t_test <- c(t_test, a)
  }
  t_test <- p.adjust(t_test)
  data_stats_amp_type$pvalue <- t_test
  data_stats_amp_type$annot <- ifelse(data_stats_amp_type$pvalue < 0.05 & data_stats_amp_type$pvalue > 0.01, "*",
                                 ifelse(data_stats_amp_type$pvalue < 0.01 & data_stats_amp_type$pvalue > 0.001, "**",
                                        ifelse(data_stats_amp_type$pvalue < 0.001, "***", "")))
  
  ### Combine stats and save
  data_stats <- rbind(data_stats_mid, data_stats_amp)
  data_stats$feature <- c("mid", "mid", "amp", "amp")
  
  data_stats_type <- rbind(data_stats_mid_type, data_stats_amp_type)
  data_stats_type$feature <- c(rep("mid", 7), rep("amp", 7))
  
  write.table(data_stats, file.path(outdir, paste0("fragment_nucleosome_features_", site, ".txt")), sep = "\t", row.names = FALSE)
  write.table(data_stats_type, file.path(outdir, paste0("fragment_nucleosome_features_", site, "_type.txt")), sep = "\t", row.names = FALSE)
  
  ### Plot features
  theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
                 axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 legend.position = "none",
                 legend.key = element_rect(fill = "white"),
                 legend.text = element_text(size = 12),
                 axis.text = element_text(size = 13),
                 axis.title = element_text(size = 13),
                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ### plotting limits
  bumper_mid <- (max(data_scores$mid, na.rm = TRUE) - min(data_scores$mid, na.rm = TRUE))/6
  bumper_amp <- (max(data_scores$amp, na.rm = TRUE) - min(data_scores$amp, na.rm = TRUE))/6
  
  ### Midpoint coverage
  data_scores <- data_scores %>%
    group_by(diag) %>%
    arrange(mid)
  fig_mid <- ggplot(data_scores, aes(diag, mid)) +
    geom_point(aes(group = diag, color = cancer_status), alpha = 0.5,
               position = position_dodge2(.75), size = 1.5, pch = 16) +
    geom_text(data = data_stats_mid, aes(x = diag, y = max(data_scores$mid, na.rm = TRUE) + bumper_mid*0.75, label = annot), size = 5) +
    xlab("") + 
    ylab("Midpoint Coverage") +
    ggtitle("Midpoint") + 
    scale_color_manual(values = c("black", "black", "red")) +
    theme +
    stat_summary(fun = match.fun(median), geom="crossbar", size = 0.25, width = 0.5, group = 1, show.legend = FALSE, color = "#FF0000", na.rm = TRUE) +
    scale_y_continuous(limits = c(min(data_scores$mid, na.rm = TRUE) - bumper_mid, max(data_scores$mid, na.rm = TRUE) + bumper_mid), expand = c(0,0))
  fig_mid
  
  ### Amplitude
  data_scores <- data_scores %>%
    group_by(diag) %>%
    arrange(amp)
  fig_amp <- ggplot(data_scores, aes(diag, amp, fill = diag)) +
    geom_point(aes(group = diag, color = cancer_status), alpha = 0.5,
               position = position_dodge2(.75), size = 1.5, pch = 16) +
    geom_text(data = data_stats_amp, aes(x = diag, y = max(data_scores$amp, na.rm = TRUE) + bumper_amp*0.75, label = annot), size = 5) +
    xlab("") + 
    ylab("Amplitude") +
    labs(color = "Cancer Status") +
    ggtitle("Amplitude") + 
    scale_color_manual(values = c("black", "black", "red")) +
    theme +
    stat_summary(fun = match.fun(median), geom="crossbar", size = 0.25, width = 0.5, group = 1, show.legend = FALSE, color = "#FF0000", na.rm = TRUE) +
    scale_y_continuous(limits = c(min(data_scores$amp, na.rm = TRUE) - bumper_amp, max(data_scores$amp, na.rm = TRUE) + bumper_amp), expand = c(0,0))
  fig_amp
  
  Figure <- ggarrange(fig_mid, fig_amp, nrow = 1)
  Figure <- annotate_figure(Figure, top = text_grob(name, size = 15, face = "bold"))
  Figure
  ggsave(file.path(outdir, paste0("fragment_nucleosome_features_", site, ".pdf")), Figure, device = "pdf", width = 3, height = 3, units = "in")
  
  ### Midpoint coverage
  data_scores <- data_scores %>%
    group_by(mutation_type) %>%
    arrange(mid)
  fig_mid_type <- ggplot(data_scores, aes(mutation_type, mid)) +
    geom_point(aes(group = diag, color = cancer_status), alpha = 0.5,
               position = position_dodge2(.75), size = 1.5, pch = 16) +
    geom_text(data = data_stats_mid_type, aes(x = mutation_type, y = min(data_scores$mid, na.rm = TRUE) - bumper_mid*0.5, label = N), size = 4) +
    geom_text(data = data_stats_mid_type, aes(x = mutation_type, y = max(data_scores$mid, na.rm = TRUE) + bumper_mid*0.75, label = annot), size = 5) +
    geom_hline(yintercept = data_stats_mid_type$Median[data_stats_mid_type$mutation_type == "Healthy"], linetype = "dashed") +
    xlab("") + 
    ylab("Midpoint Coverage") +
    ggtitle("Midpoint") + 
    scale_color_manual(values = c("black", "black", "red")) +
    theme +
    stat_summary(fun = match.fun(median), geom="crossbar", size = 0.25, width = 0.5, group = 1, show.legend = FALSE, color = "#FF0000", na.rm = TRUE) +
    scale_y_continuous(limits = c(min(data_scores$mid, na.rm = TRUE) - bumper_mid, max(data_scores$mid, na.rm = TRUE) + bumper_mid), expand = c(0,0))
  fig_mid_type
  
  ### Amplitude
  data_scores <- data_scores %>%
    group_by(mutation_type) %>%
    arrange(amp)
  fig_amp_type <- ggplot(data_scores, aes(mutation_type, amp)) +
    geom_point(aes(group = diag, color = cancer_status), alpha = 0.5,
               position = position_dodge2(.75), size = 1.5, pch = 16) +
    geom_text(data = data_stats_amp_type, aes(x = mutation_type, y = min(data_scores$amp, na.rm = TRUE) - bumper_amp*0.5, label = N), size = 4) +
    geom_text(data = data_stats_amp_type, aes(x = mutation_type, y = max(data_scores$amp, na.rm = TRUE) + bumper_amp*0.75, label = annot), size = 5) +
    geom_hline(yintercept = data_stats_amp_type$Median[data_stats_amp_type$mutation_type == "Healthy"], linetype = "dashed") +
    xlab("") + 
    ylab("Amplitude") +
    ggtitle("Amplitude") + 
    scale_color_manual(values = c("black", "black", "red")) +
    theme +
    stat_summary(fun = match.fun(median), geom="crossbar", size = 0.25, width = 0.5, group = 1, show.legend = FALSE, color = "#FF0000", na.rm = TRUE) +
    scale_y_continuous(limits = c(min(data_scores$amp, na.rm = TRUE) - bumper_amp, max(data_scores$amp, na.rm = TRUE) + bumper_amp), expand = c(0,0))
  fig_amp_type
  
  Figure <- ggarrange(fig_mid_type, fig_amp_type, nrow = 1)
  Figure <- annotate_figure(Figure, top = text_grob(name, size = 15, face = "bold"))
  Figure
  ggsave(file.path(outdir, paste0("fragment_nucleosome_features_", site, "_type.pdf")), Figure, device = "pdf", width = 6, height = 3, units = "in")
}







