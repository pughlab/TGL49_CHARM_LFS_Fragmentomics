library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpmisc)
library(ggprism)
library(ggpubr)
library(GGally)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/peak_reads"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

sites <- c("DHS", "Immune")
names <- c("Closed-Associated\nChromatin", "Open-Associated\nChromatin")

data_medians <- c()
data_proportions <- c()
data_fractions <- c()

for (i in c(1:2)) {
  ### Set variables
  site <- sites[[i]]
  name <- names[[i]]
  
  ### Find files
  data_site <- read.delim(list.files(file.path(path, "peak_reads"), paste0("freq_", site), full.names = TRUE))
  data_site_counts <- read.delim(list.files(file.path(path, "peak_reads"), paste0("fragments_", site), full.names = TRUE))
  data_freq <- read.delim(list.files(file.path(path, "insert_size"), "LFS_fragment_freq", full.names = TRUE))
  data_freq_counts <- read.delim(list.files(file.path(path, "insert_size"), "LFS_fragment.txt", full.names = TRUE))
  normal_site <- read.delim(list.files(file.path(healthy_path, "peak_reads"), paste0("freq_", site), full.names = TRUE))
  normal_site_counts <- read.delim(list.files(file.path(healthy_path, "peak_reads"), paste0("fragments_", site), full.names = TRUE))
  normal_freq <- read.delim(list.files(file.path(healthy_path, "insert_size"), "freq", full.names = TRUE))
  normal_freq_counts <- read.delim(list.files(file.path(healthy_path, "insert_size"), "fragment.txt", full.names = TRUE))
  data_samples <- read.delim(samples)
  
  ### Remove failed and unknown samples and format 
  exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0209_Cf_U_PE_373_WG")
  abnormal <- c()
  data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude, abnormal)), ]
  data_samples <- data_samples[data_samples$sWGS %in% colnames(data_site), ]
  data_freq <- data_freq[, !(colnames(data_freq) %in% c(exclude, abnormal))]
  data_freq_counts <- data_freq_counts[, !(colnames(data_freq_counts) %in% c(exclude, abnormal))]
  data_site <- data_site[, !(colnames(data_site) %in% c(exclude, abnormal))]
  data_site_counts <- data_site_counts[, !(colnames(data_site_counts) %in% c(exclude, abnormal))]
  
  ### Order samples
  common_samples <- data_samples$sWGS[data_samples$sWGS %in% colnames(data_site)]
  common_samples <- common_samples[common_samples %in% colnames(data_freq)]
  data_site <- data_site[, c("length", common_samples)]
  data_site_counts <- data_site_counts[, c("length", common_samples)]
  data_freq <- data_freq[, c("length", common_samples)]
  data_freq_counts <- data_freq_counts[, c("length", common_samples)]
  
  common_samples <- colnames(normal_site)[colnames(normal_site) %in% colnames(normal_freq)]
  normal_site <- normal_site[, common_samples]
  normal_site_counts <- normal_site_counts[, common_samples]
  normal_freq <- normal_freq[, common_samples]
  normal_freq_counts <- normal_freq_counts[, common_samples]
  
  ### Make median tables
  lfs_median <- data.frame(length = data_freq$length,
                           median = rowMedians(as.matrix(data_freq[,2:ncol(data_freq)])),
                           sd = rowSds(as.matrix(data_freq[,2:ncol(data_freq)])),
                           diag = "LFS", type = "Genome")
  lfs_site_median <- data.frame(length = data_site$length,
                                median = rowMedians(as.matrix(data_site[,2:ncol(data_site)])),
                                sd = rowSds(as.matrix(data_site[,2:ncol(data_site)])),
                                diag = "LFS", type = name)
  hbc_median <- data.frame(length = normal_freq$length,
                           median = rowMedians(as.matrix(normal_freq[,2:ncol(normal_freq)])),
                           sd = rowSds(as.matrix(normal_freq[,2:ncol(normal_freq)])),
                           diag = "HBC", type = "Genome")
  hbc_site_median <- data.frame(length = normal_site$length,
                                median = rowMedians(as.matrix(normal_site[,2:ncol(normal_site)])),
                                sd = rowSds(as.matrix(normal_site[,2:ncol(normal_site)])),
                                diag = "HBC", type = name)
  medians <-bind_rows(lfs_median, lfs_site_median, hbc_median, hbc_site_median)
  data_medians <- rbind(data_medians, medians)
  
  ### Calculate proportions
  short <- c(1:89) 
  
  lfs_prop <- data.frame(sample = colnames(data_freq[,2:ncol(data_freq)]),
                         short = colSums(data_freq[data_freq$length %in% short, 2:ncol(data_freq)])/colSums(data_freq[, 2:ncol(data_freq)]),
                         diag = "LFS", type = "Genome")
  lfs_site_prop <- data.frame(sample = colnames(data_site[,2:ncol(data_site)]),
                              short = colSums(data_site[data_site$length %in% short, 2:ncol(data_site)])/colSums(data_site[, 2:ncol(data_site)]),
                              diag = "LFS", type = name)
  hbc_prop <- data.frame(sample = colnames(normal_freq[,2:ncol(normal_freq)]),
                         short = colSums(normal_freq[normal_freq$length %in% short, 2:ncol(normal_freq)])/colSums(normal_freq[, 2:ncol(normal_freq)]),
                         diag = "HBC", type = "Genome")
  hbc_site_prop <- data.frame(sample = colnames(normal_site[,2:ncol(normal_site)]),
                              short = colSums(normal_site[normal_site$length %in% short, 2:ncol(normal_site)])/colSums(normal_site[, 2:ncol(normal_site)]),
                              diag = "HBC", type = name)
  
  proportions <- bind_rows(lfs_prop, lfs_site_prop, hbc_prop, hbc_site_prop)
  proportions$type <- factor(proportions$type, levels = c("Genome", name))
  data_proportions <- rbind(data_proportions, proportions)
  
  ### Calculate fraction of reads spanning site
  data_site_counts <- colSums(as.matrix(data_site_counts[, -1]))
  data_freq_counts <- colSums(as.matrix(data_freq_counts[,-1]))
  normal_site_counts <- colSums(as.matrix(normal_site_counts[,-1]))
  normal_freq_counts <- colSums(as.matrix(normal_freq_counts[,-1]))
  
  lfs_fraction <- data.frame(sample = colnames(data_site[,2:ncol(data_site)]),
                             fraction = data_site_counts/data_freq_counts*1000,
                             diag = "LFS", type = name)
  normal_fraction <- data.frame(sample = colnames(normal_site[,2:ncol(normal_site)]),
                                fraction = normal_site_counts/normal_freq_counts*1000,
                                diag = "HBC", type = name)
  fractions <- bind_rows(lfs_fraction, normal_fraction)
  fractions$type <- factor(fractions$type, levels = c("Genome", name))
  data_fractions <- rbind(data_fractions, fractions)
}

### Format plotting dataframes
data_medians$diag <- factor(data_medians$diag, levels = c("HBC", "LFS"),
                            labels = c("Healthy", "LFS"))
data_medians$type <- factor(data_medians$type, levels = c("Genome", "Closed-Associated\nChromatin", "Open-Associated\nChromatin"))

data_proportions$diag <- factor(data_proportions$diag, levels = c("HBC", "LFS"),
                            labels = c("Healthy", "LFS"))
data_proportions$type <- factor(data_proportions$type, levels = c("Genome", "Closed-Associated\nChromatin", "Open-Associated\nChromatin"))

data_fractions$diag <- factor(data_fractions$diag, levels = c("HBC", "LFS"),
                              labels = c("Healthy", "LFS"))

### Calculate statistics
### (proportions)
data_stats_prop1 <- data_proportions %>%
  rstatix::group_by(type) %>%
  rstatix::t_test(short~diag) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "type")
data_stats_prop1$p.adj.signif <- ifelse(data_stats_prop1$p.adj.signif == "ns", NA, data_stats_prop1$p.adj.signif)
data_stats_prop1$xmin <- c(0.8, 1, 1.7)
data_stats_prop1$xmax <- c(1.8, 2, 2.7)
data_stats_prop1$yposition <- c(0.082, 0.089, NA)

data_stats_prop2 <- data_proportions %>%
  rstatix::group_by(diag) %>%
  rstatix::t_test(short~type) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "diag")
data_stats_prop2 <- data_stats_prop2[data_stats_prop2$group1 == "Genome", ]
data_stats_prop2$p.adj.signif <- ifelse(data_stats_prop2$p.adj.signif == "****", "***", data_stats_prop2$p.adj.signif)
data_stats_prop2$p.adj.signif <- ifelse(data_stats_prop2$p.adj.signif == "ns", NA, data_stats_prop2$p.adj.signif)
data_stats_prop2$yposition <- c(NA, NA, NA, 0.075)

### (Fractions)
data_stats_frac <- data_fractions %>%
  group_by(type, diag) %>%
  dplyr::summarise(mean = mean(fraction),
                   median = median(fraction),
                   sd = sd(fraction))
a <- t.test(data_fractions$fraction[data_fractions$diag == "Healthy" & data_fractions$type == "Closed-Associated\nChromatin"], 
            data_fractions$fraction[data_fractions$diag == "LFS" & data_fractions$type == "Closed-Associated\nChromatin"])$p.value
b <- t.test(data_fractions$fraction[data_fractions$diag == "Healthy" & data_fractions$type == "Open-Associated\nChromatin"], 
            data_fractions$fraction[data_fractions$diag == "LFS" & data_fractions$type == "Open-Associated\nChromatin"])$p.value
data_stats_frac$t_test <- c(1,a,1,b)
data_stats_frac$annot <- ifelse(data_stats_frac$t_test < 0.05 & data_stats_frac$t_test > 0.01, "*",
                                ifelse(data_stats_frac$t_test < 0.01 & data_stats_frac$t_test > 0.001, "**",
                                       ifelse(data_stats_frac$t_test < 0.001, "***", NA)))

### Plot data
### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "right",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 12),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               strip.background = element_blank(),
               strip.text = element_text(size = 12))

### Frequency distributions
fig_freq <- ggplot(data_medians) +
  geom_line(aes(length, median, color = diag), alpha = 0.5) +
  geom_ribbon(aes(length, ymin = median - sd, ymax = median + sd, fill = diag), alpha = 0.1) +
  geom_vline(xintercept = 90, linetype = "dashed", size = 0.25) +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red"),
                    guide = "none") +
  facet_grid(.~type) +
  labs(fill = "", color = "Region") +
  xlab("Fragment Length (bp)") + 
  ylab("Frequency (%)") +
  ggtitle("") + 
  theme +
  theme(legend.position = c(0.2, 0.7)) +
  scale_x_continuous(limits = c(35, 150)) +
  scale_y_continuous(limits = c(-0.1, 1))
fig_freq 

ggsave(file.path(outdir, paste0("peak_reads_distributions.pdf")), height = 4, width = 8)

fig_prop <- ggplot(data_proportions) +
  geom_boxplot(aes(diag, short, fill = type), outlier.shape = 20) +
  scale_fill_manual(labels = c("Genome", "Closed", "Open"),
                    values = c("grey65", "#FB9A99", "#A6CEE3")) +
  labs(fill = "") +
  xlab("") + 
  ylab("Proportion under 90bp") +
  ggtitle("Proportion of Short Fragments") + 
  theme +
  add_pvalue(data_stats_prop1, xmin = "xmin", xmax = "xmax", y.position = "yposition", label = "p.adj.signif", tip.length = 0, label.size = 6) +
  add_pvalue(data_stats_prop2, xmin = "xmin", xmax = "xmax", y.position = "yposition", label = "p.adj.signif", tip.length = 0, label.size = 6)
fig_prop

data_fractions$type <- factor(data_fractions$type, levels = c("Closed-Associated\nChromatin", "Open-Associated\nChromatin"),
                              labels = c("Closed", "Open"))
data_stats_frac$type <- factor(data_stats_frac$type, levels = c("Closed-Associated\nChromatin", "Open-Associated\nChromatin"),
                               labels = c("Closed", "Open"))
fig_frac <- ggplot(data_fractions) +
  geom_boxplot(aes(diag, fraction, fill = diag), outlier.shape = 20) +
  geom_text(data = data_stats_frac, aes(x = diag, y = 3.5, label = annot), size = 6) +
  scale_fill_manual(values = c("grey65", "#FB9A99")) +
  facet_grid(.~type) +
  labs(fill = "") +
  xlab("") + 
  ylab("Normalized Read Fraction") +
  ggtitle("Fraction of Reads Spanning Site") + 
  theme
fig_frac

figure <- ggarrange(fig_prop, fig_frac, nrow = 1, widths = c(1, 0.8), align = "h")
figure
ggsave(file.path(outdir, paste0("peak_reads.pdf")), height = 4, width = 8)







