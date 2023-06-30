library(tidyverse)
library(dplyr)
library(matrixStats)
library(cowplot)
library(ggh4x)
library(zplyr)
library(patchwork)
library(ggpubr)
library(chisq.posthoc.test)

### Set paths
path <- "data/end_motifs"
outdir <- ""
healthy_path <- "hbc/end_motifs"

### Find paths
data <- list.files(path, "motifs.txt", full.names = TRUE)
data_lo <- data[grepl("DNASE1L3", data)]
data <- data[grepl("genome2", data)]
data_normal <- list.files(healthy_path, "motifs.txt", full.names = TRUE)
data_normal <- data_normal[grepl("genome2", data_normal)]

### Import data 
data <- read.delim(data)
data_normal <- read.delim(data_normal)
data_samples <- read.delim("sample_list.txt")

data_lo <- read.delim(data_lo)
samples_lo <- read.delim("DNASE1L3_list.txt")

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
#abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]

### Combine data
data_melt <- reshape2::melt(data, id = "motif")
data_melt <- merge(data_melt, data_samples, by.x = "variable", by.y = "sWGS")
data_melt$diag <- ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "no", "LFS-H",
                         ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "yes", "LFS-PC",
                                ifelse(data_melt$cancer_status == "positive", "LFS-AC", "Healthy")))

normal_melt <- reshape2::melt(data_normal, id = "motif")
normal_melt$diag <- "Healthy"

data_melt <- bind_rows(data_melt, normal_melt)
data_melt <- data_melt[, c("variable", "motif", "value", "diag")]

data_melt$value <- data_melt$value*100
data_melt$diag <- factor(data_melt$diag, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))

### Format DNASE1L3 data
data_lo_melt <- reshape2::melt(data_lo, id = "motif")
data_lo_melt <- merge(data_lo_melt, samples_lo, by.x = "variable", by.y = "sample")
data_lo_melt <- data_lo_melt[!(data_lo_melt$group == "SLE"), ]
data_lo_melt$group <- factor(data_lo_melt$group, levels = c("healthy", "1L3_relative", "positive"),
                             labels = c("HBC", "Het", "Homo"))

### Calculate stats
data_stats <- data_melt %>%
  group_by(motif, diag) %>%
  dplyr::summarise(mean=mean(value),
                   sd=sd(value))

motifs <- unique(data_stats$motif)
patients <- unique(data_stats$diag)
x <- c()
y <- c()
for (motif in motifs) {
  for (patient in patients) {
    a <- t.test(data_melt$value[data_melt$diag == "Healthy" & data_melt$motif == motif], 
                data_melt$value[data_melt$diag == patient & data_melt$motif == motif])$p.value
    x <- c(x, a)
    
    b <- mean(data_melt$value[data_melt$diag == patient & data_melt$motif == motif])/mean(data_melt$value[data_melt$diag == "Healthy" & data_melt$motif == motif])
    y <- c(y, b)
  }
}
x <- p.adjust(x)
data_stats$pvalue <- x
data_stats$foldchange <- y

### Calculate change of DNASE1L3 motifs
data_dnase <- data_melt[data_melt$motif %in% c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT"), ]

data_stats_dnase <- data_dnase %>%
  group_by(motif, diag) %>%
  dplyr::summarise(mean=mean(value),
                   sd=sd(value))

motifs <- unique(data_stats_dnase$motif)
patients <- unique(data_stats_dnase$diag)
x <- c()
y <- c()
for (motif in motifs) {
  for (patient in patients) {
    a <- t.test(data_dnase$value[data_dnase$diag == "Healthy" & data_dnase$motif == motif], 
                data_dnase$value[data_dnase$diag == patient & data_dnase$motif == motif])$p.value
    x <- c(x, a)
    
    b <- mean(data_dnase$value[data_dnase$diag == patient & data_dnase$motif == motif])/mean(data_dnase$value[data_dnase$diag == "Healthy" & data_dnase$motif == motif])
    y <- c(y, b)
  }
}
data_stats_dnase$pvalue <- x
data_stats_dnase$foldchange <- y
data_stats_dnase$annot <- ifelse(data_stats_dnase$pvalue < 0.05 & data_stats_dnase$pvalue > 0.01, "*",
                             ifelse(data_stats_dnase$pvalue < 0.01 & data_stats_dnase$pvalue > 0.001, "**",
                                    ifelse(data_stats_dnase$pvalue < 0.001, "***", "")))

### Calculate change of DNASE1L3 motifs (DNASE1L3 data)
data_dnase_lo <- data_lo_melt[data_lo_melt$motif %in% c("CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT"), ]

data_stats_lo <- data_dnase_lo %>%
  group_by(motif, group) %>%
  dplyr::summarise(mean=mean(value),
                   sd=sd(value))

motifs <- unique(data_stats_lo$motif)
patients <- unique(data_stats_lo$group)
x <- c()
y <- c()
for (motif in motifs) {
  for (patient in patients) {
    a <- t.test(data_dnase_lo$value[data_dnase_lo$group == "HBC" & data_dnase_lo$motif == motif], 
                data_dnase_lo$value[data_dnase_lo$group == patient & data_dnase_lo$motif == motif])$p.value
    x <- c(x, a)
    
    b <- mean(data_dnase_lo$value[data_dnase_lo$group == patient & data_dnase_lo$motif == motif])/mean(data_dnase_lo$value[data_dnase_lo$group == "HBC" & data_dnase_lo$motif == motif])
    y <- c(y, b)
  }
}
data_stats_lo$pvalue <- x
data_stats_lo$foldchange <- y
data_stats_lo$annot <- ifelse(data_stats_lo$pvalue < 0.05 & data_stats_lo$pvalue > 0.01, "*",
                                 ifelse(data_stats_lo$pvalue < 0.01 & data_stats_lo$pvalue > 0.001, "**",
                                        ifelse(data_stats_lo$pvalue < 0.001, "***", "")))

### Calculate median fold changes
motifs <- data$motif
median_hbc <- rowMedians(as.matrix(data_normal[,-1]))
median_neg <- rowMedians(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no"])]))
median_sur <- rowMedians(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes"])]))
median_pos <- rowMedians(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "positive"])]))

median_lfs_neg <- rowMedians(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "negative"])]))
sd_lfs_neg <- rowSds(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "negative"])]))

fold_neg <- median_neg/median_hbc
fold_sur <- median_sur/median_hbc
fold_pos <- median_pos/median_hbc
fold_pos2 <- median_pos/median_neg
fold_lfs <- median_lfs_neg/median_hbc
sd_lfs_neg <- sd_lfs_neg/median_hbc

fold_change <- data.frame(motif = motifs,
                          fold_neg = fold_lfs,
                          sd = sd_lfs_neg)
fold_change <- fold_change[order(fold_change$fold_neg, decreasing = TRUE), ]
order <- fold_change$motif
fold_change$motif <- factor(fold_change$motif, levels = order)

### Make motif table
data_motif <- as.data.frame(str_split_fixed(order, "", 4))
colnames(data_motif) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4")
data_motif$bin <- c(1:nrow(data_motif))
data_motif[data_motif == "A" | data_motif == "T"] <- "A/T"
data_motif[data_motif == "G" | data_motif == "C"] <- "G/C"
data_motif <- reshape2::melt(data_motif, id = "bin")

### Make frequency change comparison table
data_freq <- data.frame(motif = motifs,
                        previvor = fold_neg,
                        cancer = fold_pos2)
data_freq$change <- ifelse(data_freq$previvor > 1, "Increase", "Decrease")
data_freq$change2 <- ifelse(data_freq$cancer > 1, "Increase", "Decrease")
data_freq$change <- factor(data_freq$change, levels = c("Increase", "Decrease"))
data_freq$change2 <- factor(data_freq$change2, levels = c("Increase", "Decrease"))
data_freq <- data_freq[order(data_freq$change,
                             data_freq$cancer, decreasing = TRUE), ]
order <- data_freq$motif
data_freq$motif <- factor(data_freq$motif, levels = order)

### Make motif table 2
data_motif2 <- as.data.frame(str_split_fixed(order, "", 4))
data_motif2$sum <- data_freq$change
colnames(data_motif2) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4", "sum")
data_motif2$bin <- c(1:nrow(data_motif2))
data_motif2[data_motif2 == "A" | data_motif2 == "T"] <- "A/T"
data_motif2[data_motif2 == "G" | data_motif2 == "C"] <- "G/C"
data_motif2 <- reshape2::melt(data_motif2, id = c("bin", "sum"))

### Make regression table
data_reg <- rbind(data.frame(var1 = fold_neg,
                             var2 = fold_sur,
                             comp = "LFS-H vs\nLFS-PC"),
                  data.frame(var1 = fold_neg,
                             var2 = fold_pos,
                             comp = "LFS-H vs\nLFS-AC"),
                  data.frame(var1 = fold_sur,
                             var2 = fold_pos,
                             comp = "LFS-PC vs\nLFS-AC"))
data_reg$comp <- factor(data_reg$comp, levels = c("LFS-H vs\nLFS-PC", "LFS-H vs\nLFS-AC", "LFS-PC vs\nLFS-AC"))

### Compare frequencies based on expected frequencies (shannon)
median_hbc <- data.frame(motif = data_normal$motif,
                         median = median_hbc)
expected_freq <- 1/256
median_hbc$expected <- ifelse(median_hbc$median > expected_freq, "Above", "Below")

data_expected_freq <- merge(median_hbc, fold_change, by = "motif")
data_expected_freq$change <- ifelse(data_expected_freq$fold_neg > 1, "Increase", "Decrease")

stat_shannon <- as.data.frame(table(data_expected_freq$expected, data_expected_freq$change))

chisq_shannon <- chisq.test(data_expected_freq$expected, data_expected_freq$change)$p.value
post_hoc <- as.data.frame(table(data_expected_freq$expected, data_expected_freq$change))
chisq_shannon_decrease <- chisq.test(post_hoc$Freq[post_hoc$Var2 == "Decrease"])$p.value
chisq_shannon_increase <- chisq.test(post_hoc$Freq[post_hoc$Var2 == "Increase"])$p.value

my_chars <- strsplit(data_expected_freq$motif, "")
char_frequency <- lapply(my_chars, function(x) {table(x)})
char_frequency <- do.call(bind_rows, char_frequency)
char_frequency <- as.data.frame(char_frequency)
char_frequency[is.na(char_frequency)] <- 0
data_expected_freq$AT <- char_frequency$A + char_frequency$T

data_expected_freq$expected <- as.factor(data_expected_freq$expected)
data_expected_freq$change <- as.factor(data_expected_freq$change)
data_expected_freq$AT <- as.factor(data_expected_freq$AT)

stat_shannon2 <- as.data.frame(table(data_expected_freq$expected, data_expected_freq$change, data_expected_freq$AT))
input <- data_expected_freq[data_expected_freq$expected == "Above", ]
input <- table(input$change, input$AT)
chisq_above <- chisq.test(input)$p.value

input <- data_expected_freq[data_expected_freq$expected == "Below", ]
input <- table(input$change, input$AT)
chisq_below <- chisq.test(input)$p.value

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

### Plot comparisons
data_stats_dnase$max <- c(rep(1.13, 4),
                          rep(1.6, 4),
                          rep(0.85, 4),
                          rep(1.6, 4),
                          rep(1.05, 4),
                          rep(1.30, 4))
fig <- ggplot(data_dnase) +
  geom_boxplot(aes(diag, value, fill = diag), outlier.size = 1, alpha = 0.5) +
  geom_text(data = data_stats_dnase, aes(diag, max, label = annot), size = 5) +
  xlab("Frequency (%)") + 
  ylab("Frequency (%)") +
  labs(color = "", fill = "") +
  scale_fill_manual(name = " ", values =c("black", "#1F78B4", "#33A02C", "#E31A1C")) +
  facet_wrap(.~motif, scales = "free_y", nrow = 1) +
  facetted_pos_scales(y = list(scale_y_continuous(limits = c(0.5, 1.15)),
                               scale_y_continuous(limits = c(0.8, 1.6)),
                               scale_y_continuous(limits = c(0.65, 0.95)),
                               scale_y_continuous(limits = c(0.7, 1.55)),
                               scale_y_continuous(limits = c(0.65, 1.1)),
                               scale_y_continuous(limits = c(0.85, 1.4)))) +
  ggtitle(paste0("DNASE1L3 motifs")) + 
  theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
fig

ggsave(file.path(outdir, paste0("fragment_end_contexts_dnase1l3.pdf")), fig, width = 8, height = 4)

### Plot fold changes
fig_fold <- ggplot(fold_change) +
  geom_point(aes(motif, fold_neg), pch = 16, alpha = 0.5, color = "#1F78B4") +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  geom_errorbar(aes(motif, ymax = fold_neg + sd, ymin = fold_neg - sd), width = 0, alpha = 0.25) +
  xlab("") + 
  ylab("Fold Change") +
  labs(color = "", fill = "") +
  #scale_color_manual(name = " ", values =c("#1F78B4", "#33A02C", "#E31A1C")) +
  ggtitle("LFS Cancer Negative (vs HBC)") + 
  theme +
  theme(legend.position = c(0.25, 0.4),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
fig_fold

fig_motif <- ggplot(data_motif) +
  geom_raster(aes(bin, variable, fill = value)) +
  scale_fill_manual(values = c("A/T" = "#FC766AFF", "G/C" = "#5B84B1FF")) +
  xlab("Tetranucleotide") +
  ylab("Position") +
  labs(fill = "") +
  theme +
  theme(legend.position = c(0.15, 1.75),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(-1,0.5,0.5,0.5), "lines")) +
  scale_x_continuous(expand = c(0,0))
fig_motif

figure <- fig_fold/fig_motif + plot_layout(heights = c(2,1))
figure
ggsave(file.path(outdir, "fragment_end_contexts_change.pdf"), figure, width = 6.25, height = 4.25)

### Plot regressions
fig_reg <- ggplot(data_reg, aes(var1, var2)) +
  geom_point(aes(color = comp), pch = 16, alpha = 0.5) +
  stat_regline_equation(label.y = 1, aes(label = ..rr.label..)) +
  facet_grid(.~comp) +
  xlab("Fold Change (Var1)") + 
  ylab("Fold Change (Var2)") +
  labs(color = "", fill = "") +
  scale_color_manual(name = " ", values =c("#1F78B4", "#33A02C", "#E31A1C")) +
  ggtitle("Fold Change\n(vs Healthy)") + 
  theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fig_reg
ggsave(file.path(outdir, "fragment_end_contexts_change_reg.pdf"), fig_reg, width = 6, height = 3)

### Plot fold changes (LFS)
fig_fold2 <- ggplot(data_freq) +
  geom_point(aes(motif, cancer, color = change2), pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  facet_grid(~change, scales = "free_x", space = "free_x") +
  xlab("") + 
  ylab("Fold Change\n(LFS-AC vs LFS-H)") +
  labs(color = "", fill = "") +
  scale_color_manual(name = " ", values =c("#E31A1C", "#1F78B4")) +
  ggtitle("Direction of Change\n(LFS-H vs Healthy)") + 
  theme +
  theme(legend.position = "none",
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
fig_fold2

fig_motif2 <- ggplot(data_motif2, aes(bin, variable, fill = value)) +
  geom_tile() +
  facet_grid(~sum, scales = "free_x", space = "free_x", drop = TRUE) +
  scale_fill_manual(values = c("A/T" = "#FC766AFF", "G/C" = "#5B84B1FF")) +
  xlab("Tetranucleotide") +
  ylab("Position") +
  labs(fill = "") +
  theme +
  theme(legend.position = c(0.65, 1.8),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(-1,0.5,0.5,0.5), "lines")) +
  scale_x_discrete(drop = TRUE)
fig_motif2

figure <- fig_fold2/fig_motif2 + plot_layout(heights = c(2,1))
figure
ggsave(file.path(outdir, "fragment_end_contexts_change_lfs.pdf"), figure, width = 6.25, height = 4.25)

### Shannon plots
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "bottom",
               legend.key = element_rect(fill = "white"),
               legend.title = element_text(size = 12),
               legend.text = element_text(size = 12),
               strip.background = element_blank(),
               strip.text = element_text(size = 13),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

plot_shan <- ggplot(stat_shannon) +
  geom_bar(aes(Var2, Freq, fill = Var2), color = NA, stat = "identity") +
  facet_grid(.~Var1, scales = "free_x", space = "free_x") +
  xlab("Change vs TP53 Wildtype") + 
  ylab("# of Motifs") +
  ggtitle("Expected Motif Frequecy") + 
  scale_fill_manual(name = " ", values = c("#5B84B1FF", "#FC766AFF")) +
  theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot_shan
ggsave(file.path(outdir, "fragment_end_shannon_expected.pdf"), plot_shan, width = 3, height = 4)

plot_shan2 <- ggplot(stat_shannon2) +
  geom_bar(aes(Var3, Freq, fill = Var2), color = NA, stat = "identity") +
  facet_grid(Var2~Var1, scales = "free_x", space = "free_x") +
  xlab("# of A/T Nucleotides") + 
  ylab("# of Motifs") +
  ggtitle("Expected Motif Frequecy") + 
  scale_fill_manual(name = " ", values = c("#5B84B1FF", "#FC766AFF")) +
  theme +
  theme(legend.position = "none",
        legend.background = element_blank())
plot_shan2
ggsave(file.path(outdir, "fragment_end_shannon_AT.pdf"), plot_shan2, width = 4, height = 4)
