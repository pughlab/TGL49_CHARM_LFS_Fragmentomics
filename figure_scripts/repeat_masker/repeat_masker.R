library(dplyr)
library(matrixStats)
library(data.table)
library(ggplot2)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/repeat_masker"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/repeat_masker"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/repeat_masker"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

## Import data
data_count <- read.delim(list.files(path, "count", full.names = TRUE))
normal_count <- read.delim(list.files(healthy_path, "count", full.names = TRUE)) 

data_frequency <- read.delim(list.files(path, "freq", full.names = TRUE))
normal_frequency <- read.delim(list.files(healthy_path, "freq", full.names = TRUE))

data_prop <- read.delim(list.files(path, "prop.txt", full.names = TRUE))
normal_prop <- read.delim(list.files(healthy_path, "prop.txt", full.names = TRUE))

data_samples <- read.delim(samples)

### Restrict size 
data_frequency <- data_frequency[data_frequency$length %in% c(10:430), ]
normal_frequency <- normal_frequency[normal_frequency$length %in% c(10:430), ]

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_frequency), ]
data_frequency <- data_frequency[ , colnames(data_frequency) %in% c("length", "type", data_samples$sWGS)]
data_frequency[is.na(data_frequency)] <- 0

data_count <- data_count[ , colnames(data_count) %in% c("length", "type", data_samples$sWGS)]
data_count[is.na(data_count)] <- 0

data_prop <- data_prop[colnames(data_prop) %in% c("type", data_samples$sWGS)]

### Make healthy median
tests <- unique(normal_frequency$type)
normal_frequency$length <- factor(normal_frequency$length, levels = c(10:430))
normal_frequency$type <- factor(normal_frequency$type, levels = tests)

normal_frequency$median <- rowMedians(as.matrix(normal_frequency[, colnames(normal_frequency) %like% "TGL"]))
normal_median <- normal_frequency[, c("length", "type", "median")]

### Make previvor median
data_frequency$length <- factor(data_frequency$length, levels = c(10:430))
data_frequency$type <- factor(data_frequency$type, levels = tests)

samples_neg <- data_samples[data_samples$cancer_status == "negative" &
                              data_samples$previous_cancer == "no", ]
data_neg <- data_frequency[ , colnames(data_frequency) %in% c("length", "type", samples_neg$sWGS)]
data_neg$median <- rowMedians(as.matrix(data_neg[, colnames(data_neg) %like% "TGL"]))
data_neg <- data_neg[, c("length", "type", "median")]

### Make cancer survivor median
samples_sur <- data_samples[data_samples$cancer_status == "negative" &
                              data_samples$previous_cancer == "yes", ]
data_sur <- data_frequency[ , colnames(data_frequency) %in% c("length", "type", samples_sur$sWGS)]
data_sur$median <- rowMedians(as.matrix(data_sur[, colnames(data_sur) %like% "TGL"]))
data_sur <- data_sur[, c("length", "type", "median")]

### Make cancer positive median
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]
data_pos <- data_frequency[ , colnames(data_frequency) %in% c("length", "type", samples_pos$sWGS)]
data_pos$median <- rowMedians(as.matrix(data_pos[, colnames(data_pos) %like% "TGL"]))
data_pos <- data_pos[, c("length", "type", "median")]

### Calculate KS-tests
p_freq <- data.frame()

for (i in c(1:length(tests))) {
  test <- tests[[i]]
  
  w <- normal_median[normal_median$type == test, ]
  x <- data_neg[data_neg$type == test, ]
  y <- data_sur[data_sur$type == test, ]
  z <- data_pos[data_pos$type == test, ]
  
  a <- round(ks.test(w$median, x$median)$p.value, 4)
  b <- round(ks.test(w$median, y$median)$p.value, 4)
  c <- round(ks.test(w$median, z$median)$p.value, 4)
  
  t <- data.frame(type = test,
                  sample = c("LFS-H", "LFS-PC", "LFS-AC"),
                  pvalue = c(a, b, c))
  
  p_freq <- rbind(p_freq, t)
  
}
p_freq$padj <- p.adjust(p_freq$pvalue)

### Make plotting table
data <- data.frame(length = as.numeric(normal_median$length),
                   type = normal_median$type, 
                   hbc = normal_median$median,
                   neg = data_neg$median,
                   sur = data_sur$median,
                   pos = data_pos$median)
data$length <- as.numeric(data$length)
data$type <- factor(data$type, levels = tests)

### Plot differences
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.title = element_text(size = 13),
               axis.line = element_line(colour = "black"),
               axis.text = element_text(size = 11),
               legend.position = "none",
               legend.text = element_text(size = 12),
               legend.key = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(size = 11))

plot_freq <- ggplot(data) +
  geom_line(aes(length, hbc, color = "black"), alpha = 0.5) +
  geom_line(aes(length, neg, color = "#1F78B4"), alpha = 0.5) +
  geom_line(aes(length, sur, color = "#33A02C"), alpha = 0.5) +
  geom_line(aes(length, pos, color = "#E31A1C"), alpha = 0.5) +
  geom_vline(xintercept = c(168), linetype = "dashed", color = "black") +
  ggtitle("Fragment Frequency Distribution") +
  xlab("Fragment Length") + 
  ylab("Frequency (%)") +
  facet_wrap(.~type, scales = "free", ncol = 6) +
  theme +
  theme(legend.position = "bottom") +
  scale_color_manual(name = " ", 
                     values =c("black" = "black", "#1F78B4" = "#1F78B4", "#33A02C" = "#33A02C", "#E31A1C" = "#E31A1C"), 
                     labels = c("Healthy","LFS-H", "LFS-PC", "LFS-AC"))

ggsave(file.path(outdir, "repeats_frequency.pdf"), width = 14.4, height = 13)

### Seperate genome-wide
genome_hbc <- normal_prop[normal_prop$type == "genome", -1]
genome_hbc <- as.data.frame(t(genome_hbc))
normal_prop <- normal_prop[!(normal_prop$type == "genome"), ]
#normal_prop[, 2:ncol(normal_prop)] <- sweep(normal_prop[, 2:ncol(normal_prop)], 2, genome_hbc, "/")

genome_lfs <- data_prop[data_prop$type == "genome", -1]
genome_lfs <- as.data.frame(t(genome_lfs))
data_prop <- data_prop[!(data_prop$type == "genome"), ]
#data_prop[, 2:ncol(data_prop)] <- sweep(data_prop[, 2:ncol(data_prop)], 2, genome_lfs, "/")

### Format prop dataframes
samples_neg <- data_samples[data_samples$cancer_status == "negative", ]
prop_neg <- data_prop[, colnames(data_prop) %in% c("type", samples_neg$sWGS)]
prop_pos <- data_prop[, colnames(data_prop) %in% c("type", samples_pos$sWGS)]

### Statistics for proportions
p_prop <- data.frame()
tests <- tests[!(tests == "genome")]

for (i in c(1:length(tests))) {
  test <- tests[[i]]
  
  x <- t(normal_prop[normal_prop$type == test, -1])
  y <- t(prop_neg[prop_neg$type == test, -1])
  z <- t(prop_pos[prop_pos$type == test, -1])
  
  a <- round(t.test(x, y)$p.value, 4)
  b <- round(t.test(x, z)$p.value, 4)
  
  t <- data.frame(type = test,
                  sample = c("Healthy", "Negative", "Positive"),
                  pvalue = c(1, a, b))
  
  p_prop <- rbind(p_prop, t)
}

p_prop$annot <- ifelse(p_prop$pvalue < 0.05 & p_prop$pvalue > 0.01, "*",
                       ifelse(p_prop$pvalue < 0.01 & p_prop$pvalue > 0.001, "**",
                              ifelse(p_prop$pvalue < 0.001, "***", "")))

### Format proportions for plotting
data_prop <- reshape2::melt(data_prop, id = "type")
data_prop$sample <- ifelse(data_prop$variable %in% colnames(prop_neg), "Negative", "Positive")

normal_prop <- reshape2::melt(normal_prop, id = "type")
normal_prop$sample <- "Healthy"

data_prop <- bind_rows(data_prop, normal_prop)

### Plot proportions
plot_prop <- ggplot(data_prop) +
  geom_boxplot(aes(sample, value, fill = sample), outlier.size = 0.25, outlier.shape = 16, alpha = 0.5) +
  geom_text(data = p_prop, aes(sample, Inf, label = annot), vjust = 1.2) +
  ggtitle("Proportion of Short Fragments") +
  xlab("Repeat Element") + 
  ylab("Proportion <150bp") +
  facet_wrap(.~type, scales = "free", nrow = 4) +
  theme +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 10)) +
  scale_fill_manual(name = "Cancer Status", values = c("black", "#1F78B4", "#E31A1C"))
#plot_prop
ggsave(file.path(outdir, "repeats_normalized_short_prop.pdf"), width = 20, height = 9)

### Find total number of fragments
data_count <- reshape2::melt(data_count, id = c("length", "type"))
data_count <- data_count %>%
  group_by(type, variable) %>%
  dplyr::summarise(count = sum(value))

normal_count <- reshape2::melt(normal_count, id = c("length", "type"))
normal_count <- normal_count %>%
  group_by(type, variable) %>%
  dplyr::summarise(count = sum(value))

genome_count_lfs <- data_count[data_count$type == "genome", ]
data_count <- data_count[!(data_count$type == "genome"), ]

genome_count_hbc <- normal_count[normal_count$type == "genome", ]
normal_count <- normal_count[!(normal_count$type == "genome"), ]

### Normalize by total count
data_count <- merge(data_count, genome_count_lfs, by = "variable")
data_count$prop <- data_count$count.x/data_count$count.y
data_count$log <- log(data_count$prop)
data_count$status <- ifelse(data_count$variable %in% samples_neg$sWGS, "Negative", "Positive")

normal_count <- merge(normal_count, genome_count_hbc, by = "variable")
normal_count$prop <- normal_count$count.x/normal_count$count.y
normal_count$log <- log(normal_count$prop)
normal_count$status <- "Healthy"

data_plot <-bind_rows(data_count, normal_count)

### Format prop dataframes
count_neg <- data_count[data_count$variable %in% c("type", samples_neg$sWGS), ]
count_pos <- data_count[data_count$variable %in% c("type", samples_pos$sWGS), ]

### Statistics for counts
p_count <- data.frame()

for (i in c(1:length(tests))) {
  test <- tests[[i]]
  
  x <- normal_count[normal_count$type.x == test, ]
  y <- count_neg[count_neg$type.x == test, ]
  z <- count_pos[count_pos$type.x == test, ]
  
  a <- round(t.test(x$prop, y$prop)$p.value, 4)
  b <- round(t.test(x$prop, z$prop)$p.value, 4)
  
  t <- data.frame(type.x = test,
                  status = c("Healthy", "Negative", "Positive"),
                  pvalue = c(1, a, b))
  
  p_count <- rbind(p_count, t)
}

p_count$annot <- ifelse(p_count$pvalue < 0.05 & p_count$pvalue > 0.01, "*",
                       ifelse(p_count$pvalue < 0.01 & p_count$pvalue > 0.001, "**",
                              ifelse(p_count$pvalue < 0.001, "***", "")))

### Plot proportions
plot_count <- ggplot(data_plot) +
  geom_boxplot(aes(status, log, fill = status), outlier.size = 0.25, outlier.shape = 16, alpha = 0.5) +
  geom_text(data = p_count, aes(status, Inf, label = annot), vjust = 1.2) +
  ggtitle("Normalized Proportion of All Fragments") +
  xlab("Repeat Element") + 
  ylab("Log10(Proportion)") +
  facet_wrap(.~type.x, scales = "free", nrow = 4) +
  theme +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 10)) +
  scale_fill_manual(name = "Cancer Status", values = c("black", "#1F78B4", "#E31A1C"))
#plot_count
ggsave(file.path(outdir, "repeats_normalized_total_prop.pdf"), width = 20, height = 9)

### Caclulate integrated contribution
data_cohort_short <- bind_rows(genome_lfs, genome_hbc)

data_prop_short <- data_prop
colnames(data_prop_short) <- c("type", "sample", "short", "status")
data_prop_short <- merge(data_prop_short, data_cohort_short, by.x = "sample", by.y = "row.names")
data_prop_short$short <- data_prop_short$short/data_prop_short$`15`

data_prop_all <- data_plot[, c("variable", "type.x", "prop", "status")]
colnames(data_prop_all) <- c("sample", "type", "all", "status")

data_contribution <- merge(data_prop_short, data_prop_all, by = c("sample", "type", "status"))
data_contribution$score <- data_contribution$short*data_contribution$all
data_contribution$log <- log(data_contribution$score)

### Format contribution dataframes
contribution_neg <- data_contribution[data_contribution$sample %in% samples_neg$sWGS, ]
contribution_pos <- data_contribution[data_contribution$sample %in% samples_pos$sWGS, ]
contribution_hbc <- data_contribution[data_contribution$status == "Healthy", ]

### Statistics for contribution scores
p_contribution <- data.frame()

for (i in c(1:length(tests))) {
  test <- tests[[i]]
  
  x <- contribution_hbc[contribution_hbc$type == test, ]
  y <- contribution_neg[contribution_neg$type == test, ]
  z <- contribution_pos[contribution_pos$type == test, ]
  
  a <- round(t.test(x$score, y$score)$p.value, 4)
  b <- round(t.test(x$score, z$score)$p.value, 4)
  
  t <- data.frame(type = test,
                  status = c("Healthy", "Negative", "Positive"),
                  pvalue = c(1, a, b))
  
  p_contribution <- rbind(p_contribution, t)
}

p_contribution$annot <- ifelse(p_contribution$pvalue < 0.05 & p_contribution$pvalue > 0.01, "*",
                        ifelse(p_contribution$pvalue < 0.01 & p_contribution$pvalue > 0.001, "**",
                               ifelse(p_contribution$pvalue < 0.001, "***", "")))

### Plot proportions
plot_contribution <- ggplot(data_contribution) +
  geom_boxplot(aes(status, log, fill = status), outlier.size = 0.25, outlier.shape = 16, alpha = 0.5) +
  geom_text(data = p_contribution, aes(status, Inf, label = annot), vjust = 1.2) +
  ggtitle("Normalized Contribution of Short Fragments") +
  xlab("Repeat Element") + 
  ylab("Log10(Contribution)") +
  facet_wrap(.~type, scales = "free", nrow = 4) +
  theme +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 10)) +
  scale_fill_manual(name = "Cancer Status", values = c("black", "#1F78B4", "#E31A1C"))
#plot_contribution
ggsave(file.path(outdir, "repeats_normalized_contribution.pdf"), width = 20, height = 9)

### Make a zoomed in plot
samples_neg <- samples_neg[samples_neg$previous_cancer == "no", ]
contribution_neg <- data_contribution[data_contribution$sample %in% samples_neg$sWGS, ]
contribution_sur <- data_contribution[data_contribution$sample %in% samples_sur$sWGS, ]

p_zoom <- data.frame()

for (i in c(1:length(tests))) {
  test <- tests[[i]]
  
  w <- contribution_hbc[contribution_hbc$type == test, ]
  x <- contribution_neg[contribution_neg$type == test, ]
  y <- contribution_sur[contribution_sur$type == test, ]
  z <- contribution_pos[contribution_pos$type == test, ]
  
  a <- round(t.test(w$score, x$score)$p.value, 4)
  b <- round(t.test(w$score, y$score)$p.value, 4)
  c <- round(t.test(w$score, z$score)$p.value, 4)
  
  t <- data.frame(type = test,
                  status = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"),
                  pvalue = c(1, a, b, c))
  
  p_zoom <- rbind(p_zoom, t)
}

p_zoom$annot <- ifelse(p_zoom$pvalue < 0.05 & p_zoom$pvalue > 0.01, "*",
                       ifelse(p_zoom$pvalue < 0.01 & p_zoom$pvalue > 0.001, "**",
                              ifelse(p_zoom$pvalue < 0.001, "***", "")))

data_zoom <- data_contribution[data_contribution$type %in% c("GAATG", "LINE_L2", "LTR", "LTR_ERV1", "LTR_ERVL", "LTR_EVRK", "LTR_ERVL.MaLR", "LTR_Gypsy",
                                                             "SINE_MIR", "SINE_tRNA", "BSR_Beta"), ]
data_zoom$type <- factor(data_zoom$type, levels = c("LTR", "LTR_ERV1", "LTR_EVRK", "LTR_ERVL", "LTR_ERVL.MaLR", "LTR_Gypsy",
                                                    "GAATG", "BSR_Beta", "LINE_L2", "SINE_MIR", "SINE_tRNA"),
                         labels = c("LTR", "LTR/ERV1", "LTR/ERVK", "LTR/ERVL", "LTR/ERVL-MaLR", "LTR/Gypsy",
                                    "(GAATG)n", "BSR/Beta", "LINE-L2", "SINE/MIR", "SINE/tRNA"))
data_zoom$status <- ifelse(data_zoom$sample %in% samples_neg$sWGS, "LFS-H", data_zoom$status)
data_zoom$status <- ifelse(data_zoom$sample %in% samples_sur$sWGS, "LFS-PC", data_zoom$status)
data_zoom$status <- ifelse(data_zoom$sample %in% samples_pos$sWGS, "LFS-AC", data_zoom$status)
data_zoom$status <- factor(data_zoom$status, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))

data_zoom_p <- p_zoom[p_zoom$type %in% c("GAATG", "LINE_L2", "LTR", "LTR_ERV1", "LTR_ERVL", "LTR_EVRK", "LTR_ERVL.MaLR", "LTR_Gypsy",
                                         "SINE_MIR", "SINE_tRNA", "BSR_Beta"), ]
data_zoom_p$type <- factor(data_zoom_p$type, levels = c("LTR", "LTR_ERV1", "LTR_EVRK", "LTR_ERVL", "LTR_ERVL.MaLR", "LTR_Gypsy",
                                                        "GAATG", "BSR_Beta", "LINE_L2", "SINE_MIR", "SINE_tRNA"),
                           labels = c("LTR", "LTR/ERV1", "LTR/ERVK", "LTR/ERVL", "LTR/ERVL-MaLR", "LTR/Gypsy",
                                      "(GAATG)n", "BSR/Beta", "LINE-L2", "SINE/MIR", "SINE/tRNA"))
data_zoom_p$status <- factor(data_zoom_p$status, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))

data_zoom <- data_zoom %>%
  group_by(type, status) %>%
  arrange(log)

plot_zoom <- ggplot(data_zoom) +
  geom_boxplot(aes(status, log, fill = status), outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(status, log, group = status), color = "black", alpha = 0.25,
             position = position_dodge2(.75), size = 1, pch = 16) +
  geom_text(data = data_zoom_p, aes(status, Inf, label = annot), size = 5, vjust = 1.2) +
  ggtitle("Normalized Contribution of Short Fragments") +
  xlab("Repeat Element") + 
  ylab("Log10(Contribution)") +
  facet_wrap(.~type, scales = "free", nrow = 2) +
  theme +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  scale_fill_manual(name = "Cancer Status", values = c("black", "#1F78B4", "#33A02C", "#E31A1C"))
plot_zoom
ggsave(file.path(outdir, "repeats_zoomed.pdf"), width = 8, height = 6)









