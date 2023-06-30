library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- "data/insert_size"
outdir <- ""
healthy_path <- "hbc/insert_size"
samples <- "sample_list.txt"

## Import data
data_frequency <- read.delim(list.files(path, "LFS_fragment_freq.txt", full.names = TRUE))
data_normal <- read.delim(list.files(healthy_path, "freq.txt", full.names = TRUE))
data_samples <- read.delim(samples)

### Restrict size 
data_frequency <- data_frequency[data_frequency$length %in% c(15:500), ]
data_normal <- data_normal[data_normal$length %in% c(15:500), ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[(data_samples$cancer_status %in% c("positive", "negative")), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_frequency), ]
row.names(data_frequency) <- data_frequency$length
data_frequency <- data_frequency[ , colnames(data_frequency) %in% data_samples$sWGS]
data_frequency[is.na(data_frequency)] <- 0

## Set samples to graph and order factors
data_samples$sWGS <- as.factor(data_samples$sWGS)
data_samples$cancer_status <- factor(data_samples$cancer_status,
                                     levels = c("positive", "negative"),
                                     labels = c("Positive", "Negative"))
row.names(data_samples) <- data_samples$sWGS

### Make healthy median
row.names(data_normal) <- data_normal$length
data_normal <- data_normal[, !(colnames(data_normal) == "length")]
normal_median <- rowMedians(as.matrix(data_normal))
normal_sd <- rowSds(as.matrix(data_normal))

### Make LFS median
lfs_median <- rowMedians(as.matrix(data_frequency))
lfs_sd <- rowSds(as.matrix(data_frequency))

### Make previvor median
samples_neg <- data_samples[data_samples$cancer_status == "Negative" &
                              data_samples$previous_cancer == "no", ]
data_neg <- data_frequency[ , colnames(data_frequency) %in% samples_neg$sWGS]
previvor_median <- rowMedians(as.matrix(data_neg))
previvor_sd <- rowSds(as.matrix(data_neg))

### Make cancer survivor median
samples_neg <- data_samples[data_samples$cancer_status == "Negative" &
                              data_samples$previous_cancer == "yes", ]
data_neg <- data_frequency[ , colnames(data_frequency) %in% samples_neg$sWGS]
negative_median <- rowMedians(as.matrix(data_neg))
negative_sd <- rowSds(as.matrix(data_neg))

### Make cancer positive median
samples_pos <- data_samples[data_samples$cancer_status == "Positive", ]
data_pos <- data_frequency[ , colnames(data_frequency) %in% samples_pos$sWGS]
positive_median <- rowMedians(as.matrix(data_pos))
positive_sd <- rowSds(as.matrix(data_pos))

### Make plotting table
data <- data.frame(length = as.numeric(row.names(data_frequency)),
                   normal_median = normal_median, 
                   previvor_median = previvor_median, 
                   negative_median = negative_median, 
                   positive_median = positive_median,
                   lfs_median = lfs_median)

### Get peaks
normal_max <- data$length[data$normal_median == max(data$normal_median)]
previvor_max <- data$length[data$previvor_median == max(data$previvor_median)]
negative_max <- data$length[data$negative_median == max(data$negative_median)]
positive_max <- data$length[data$positive_median == max(data$positive_median)]

### Plot differences
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.title = element_text(size = 13),
               axis.line = element_line(colour = "black"),
               axis.text = element_text(size = 13),
               legend.position = "none",
               legend.text = element_text(size = 13),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank())

plot_freq <- ggplot(data) +
  geom_line(aes(length, normal_median, color = "black")) +
  geom_line(aes(length, previvor_median, color = "#1F78B4")) +
  geom_line(aes(length, negative_median, color = "#33A02C")) +
  geom_line(aes(length, positive_median, color = "#E31A1C")) +
  geom_vline(xintercept = c(168), linetype = "dashed", color = "black") +
  ggtitle("Fragment Frequency Distribution") +
  xlab("Fragment Length") + 
  ylab("Frequency (%)") +
  theme +
  theme(legend.position = "bottom",
        legend.key = element_blank()) +
  scale_color_manual(name = " ", 
                     values =c("black" = "black", "#1F78B4" = "#1F78B4", "#33A02C" = "#33A02C", "#E31A1C" = "#E31A1C"), 
                     labels = c("Healthy","LFS-H", "LFS-PC", "LFS-AC")) +
  scale_x_continuous(limits = c(0,400)) +
  scale_y_continuous(limits=c(0, 3.5), expand = c(0,0))
plot_freq

plot_inset <- ggplot(data) +
  geom_line(aes(length, normal_median, color = "black")) +
  geom_line(aes(length, previvor_median, color = "#1F78B4")) +
  geom_line(aes(length, negative_median, color = "#33A02C")) +
  geom_line(aes(length, positive_median, color = "#E31A1C")) +
  geom_vline(xintercept = c(145, 135, 123, 112, 102, 92, 82), linetype = "dashed", color = "black", size = 0.25) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  theme +
  theme(legend.position = "none",
        legend.key = element_blank(),
        plot.margin = margin(-0.5, 0, -0.5, -0.5, "cm"),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  scale_color_manual(name = " ", 
                     values =c("black" = "black", "#1F78B4" = "#1F78B4", "#33A02C" = "#33A02C", "#E31A1C" = "#E31A1C"), 
                     labels = c("Healthy","LFS-H", "LFS-PC", "LFS-AC")) +
  scale_x_continuous(limits = c(90, 150), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
plot_inset

plot_inset2 <- ggplot(data) +
  geom_line(aes(length, normal_median, color = "black")) +
  geom_line(aes(length, previvor_median, color = "#1F78B4")) +
  geom_line(aes(length, negative_median, color = "#33A02C")) +
  geom_line(aes(length, positive_median, color = "#E31A1C")) +
  geom_vline(xintercept = c(334), linetype = "dashed", color = "black", size = 0.25) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  theme +
  theme(legend.position = "none",
        legend.key = element_blank(),
        plot.margin = margin(-0.5, 0, -0.5, -0.5, "cm"),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  scale_color_manual(name = " ", 
                     values =c("black" = "black", "#1F78B4" = "#1F78B4", "#33A02C" = "#33A02C", "#E31A1C" = "#E31A1C"), 
                     labels = c("Healthy","LFS-H", "LFS-PC", "LFS-AC")) +
  scale_x_continuous(limits = c(250, 400), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0))
plot_inset2

plot_combined <- ggdraw() +
  draw_plot(plot_freq) +
  draw_plot(plot_inset, x = 0.075, y = 0.3, width = 0.3, height = 0.625) +
  draw_plot(plot_inset2, x = 0.6, y = 0.3, width = 0.35, height = 0.625)
plot_combined

ggsave(file.path(outdir, "fragment_frequency.pdf"), width = 8, height = 5)
