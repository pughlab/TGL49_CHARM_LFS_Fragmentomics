library(tidyverse)
library(dplyr)
library(matrixStats)
library(cowplot)
library(ggh4x)
library(zplyr)
library(patchwork)
library(ggpubr)

### Set paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/end_motifs"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/end_motif"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/end_motifs"

### Find paths
data <- list.files(path, "motifs.txt", full.names = TRUE)
data_lo <- data[grepl("DNASE1L3", data)]
data <- data[grepl("genome2", data)]
data_normal <- list.files(healthy_path, "motifs.txt", full.names = TRUE)
data_normal <- data_normal[grepl("genome2", data_normal)]

### Import data 
data <- read.delim(data)
data_normal <- read.delim(data_normal)
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

data_lo <- read.delim(data_lo)
samples_lo <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/DNASE1L3_list.txt")

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
#abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]

### Calculate median fold changes
motifs <- data$motif
median_hbc <- rowMedians(as.matrix(data_normal[,-1]))
median_neg <- rowMedians(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no"])]))
median_sur <- rowMedians(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes"])]))
median_pos <- rowMedians(as.matrix(data[, colnames(data) %in% unique(data_samples$sWGS[data_samples$cancer_status == "positive"])]))
median_lfs <- rowMedians(as.matrix(data[,-1]))

fold_neg <- median_neg/median_hbc
fold_sur <- median_sur/median_hbc
fold_pos <- median_pos/median_hbc
fold_lfs <- median_lfs/median_hbc

data_lo <- data_lo[order(factor(data_lo$motif, levels = motifs)), ]
median_con <- rowMedians(as.matrix(data_lo[ ,colnames(data_lo) %in% unique(samples_lo$sample[samples_lo$group == "healthy"])]))
median_dna <- rowMedians(as.matrix(data_lo[ ,colnames(data_lo) %in% unique(samples_lo$sample[samples_lo$group == "positive"])]))
median_het <- rowMedians(as.matrix(data_lo[ ,colnames(data_lo) %in% unique(samples_lo$sample[samples_lo$group == "1L3_relative"])]))

fold_dna <- median_dna/median_con
fold_het <- median_het/median_con

### Make a comparison table
data_comp <- rbind(data.frame(var1 = fold_neg,
                              var2 = fold_dna,
                              var3 = fold_het,
                              comp1 = "LFS",
                              comp2 = "DNASE",
                              comp3 = "HET"))
data_comp$ch1 <- ifelse(data_comp$var1 > 1, 1, 0)
data_comp$ch2 <- ifelse(data_comp$var2 > 1, 1, 0)
data_comp$ch3 <- ifelse(data_comp$var3 > 1, 1, 0)
data_comp$sum1 <- data_comp$ch1 + data_comp$ch2
data_comp$sum2 <- data_comp$ch1 + data_comp$ch3
data_comp$sum1 <- factor(data_comp$sum1, levels = c(2,1,0),
                         labels = c("Increase", "Discordant", "Decrease"))
data_comp$sum2 <- factor(data_comp$sum2, levels = c(2,1,0),
                         labels = c("Increase", "Discordant", "Decrease"))
data_comp$motif <- motifs

### Sum the changes
data_change_dna <- as.data.frame(table(data_comp$sum1))
data_change_het <- as.data.frame(table(data_comp$sum2))

### Format plotting tables and make motif tables (homozygous)
data_comp_dna <- data_comp[order(data_comp$sum1,
                                 data_comp$var2, decreasing = TRUE), ]
order <- data_comp_dna$motif
data_melt_dna <- reshape2::melt(data_comp_dna[, c("var1", "var2", "sum1", "motif")], id = c("motif", "sum1"))
data_melt_dna$variable <- factor(data_melt_dna$variable, levels = c("var1", "var2"),
                             labels = c("TP53 +/-", "DNASE1L3 -/-"))
data_melt_dna$motif <- factor(data_melt_dna$motif, levels = order)
data_melt_dna$value <- log2(data_melt_dna$value)

### Make motif table
data_motif_dna <- as.data.frame(str_split_fixed(order, "", 4))
data_motif_dna$sum <- data_comp_dna$sum1
colnames(data_motif_dna) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4", "sum")
data_motif_dna$bin <- c(1:nrow(data_motif_dna))
data_motif_dna[data_motif_dna == "A" | data_motif_dna == "T"] <- "A/T"
data_motif_dna[data_motif_dna == "G" | data_motif_dna == "C"] <- "G/C"
data_motif_dna <- reshape2::melt(data_motif_dna, id = c("bin", "sum"))

### Format plotting tables and make motif tables (heterozygous)
data_comp_het <- data_comp[order(data_comp$sum2,
                                 data_comp$var3, decreasing = TRUE), ]
order <- data_comp_het$motif
data_melt_het <- reshape2::melt(data_comp_het[, c("var1", "var3", "sum2", "motif")], id = c("motif", "sum2"))
data_melt_het$variable <- factor(data_melt_het$variable, levels = c("var1", "var3"),
                                 labels = c("TP53 +/-", "DNASE1L3 +/-"))
data_melt_het$motif <- factor(data_melt_het$motif, levels = order)
data_melt_het$value <- log2(data_melt_het$value)

### Make motif table
data_motif_het <- as.data.frame(str_split_fixed(order, "", 4))
data_motif_het$sum <- data_comp_het$sum2
colnames(data_motif_het) <- c("5' Base 1", "5' Base 2", "5' Base 3", "5' Base 4", "sum")
data_motif_het$bin <- c(1:nrow(data_motif_het))
data_motif_het[data_motif_het == "A" | data_motif_het == "T"] <- "A/T"
data_motif_het[data_motif_het == "G" | data_motif_het == "C"] <- "G/C"
data_motif_het <- reshape2::melt(data_motif_het, id = c("bin", "sum"))

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
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Plot number of changed motifs
data_change_dna$comp <- "TP53 +/- vs\nDNASE1L3 -/-"
data_change_het$comp <- "TP53 +/- vs\nDNASE1L3 +/-"
data_change <- bind_rows(data_change_dna, data_change_het)
fig_comp <- ggplot(data_change, aes(Var1, Freq, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.5) +
  geom_text(aes(label = Freq), vjust = -0.5) +
  facet_grid(.~comp) +
  xlab("Direction") + 
  ylab("# of Motifs") +
  labs(color = "", fill = "") +
  scale_fill_manual(name = " ", values =c("#E31A1C", "black", "#1F78B4")) +
  ggtitle("") + 
  theme +
  theme(plot.title = element_text(face = "italic"),
        panel.border = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,150))
fig_comp
ggsave(file.path(outdir, "fragment_end_contexts_dnase1l3_direction.pdf"), fig_comp, width = 3.5, height = 5)

### Plot fold changes
fig_fold <- ggplot(data_melt_dna) +
  geom_point(aes(motif, value, color = variable), pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  facet_grid(~sum1, scales = "free_x", space = "free_x") +
  xlab("") + 
  ylab("Log2(Fold Change)") +
  labs(color = "", fill = "") +
  scale_color_manual(name = " ", values =c("#1F78B4", "#33A02C", "#E31A1C")) +
  ggtitle("5' Fragment End Motifs") + 
  theme +
  theme(legend.position = c(0.55, 0.85),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.margin = margin(c(-1,2,1,1)),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
fig_fold

fig_motif <- ggplot(data_motif_dna, aes(bin, variable, fill = value)) +
  geom_tile() +
  facet_grid(~sum, scales = "free_x", space = "free_x", drop = TRUE) +
  scale_fill_manual(values = c("A/T" = "#FC766AFF", "G/C" = "#5B84B1FF")) +
  xlab("Tetranucleotide") +
  ylab("Position") +
  labs(fill = "") +
  theme +
  theme(legend.position = c(0.9, 2.9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(-1,0.5,0.5,0.5), "lines")) +
  scale_x_discrete(drop = TRUE)
fig_motif

figure <- fig_fold/fig_motif + plot_layout(heights = c(2,1))
figure
ggsave(file.path(outdir, "fragment_end_contexts_dnase1l3_change_hom.pdf"), figure, width = 7.5, height = 5)

### Plot fold changes (het)
fig_fold <- ggplot(data_melt_het) +
  geom_point(aes(motif, value, color = variable), pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  facet_grid(~sum2, scales = "free_x", space = "free_x") +
  xlab("") + 
  ylab("Log2(Fold Change)") +
  labs(color = "", fill = "") +
  scale_color_manual(name = " ", values =c("#1F78B4", "#33A02C", "#E31A1C")) +
  ggtitle("5' Fragment End Motifs") + 
  theme +
  theme(legend.position = c(0.52, 0.3),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.margin = margin(c(-1,2,1,1)),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
fig_fold

fig_motif <- ggplot(data_motif_het, aes(bin, variable, fill = value)) +
  geom_tile() +
  facet_grid(~sum, scales = "free_x", space = "free_x", drop = TRUE) +
  scale_fill_manual(values = c("A/T" = "#FC766AFF", "G/C" = "#5B84B1FF")) +
  xlab("Tetranucleotide") +
  ylab("Position") +
  labs(fill = "") +
  theme +
  theme(legend.position = c(0.9, 1.75),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(-1,0.5,0.5,0.5), "lines")) +
  scale_x_discrete(drop = TRUE)
fig_motif

figure <- fig_fold/fig_motif + plot_layout(heights = c(2,1))
figure
ggsave(file.path(outdir, "fragment_end_contexts_dnase1l3_change_het.pdf"), figure, width = 7.5, height = 5)

### Plot regressions
data_cast_dna <- reshape2::dcast(data_melt_dna, motif ~ variable)
fig_reg1 <- ggplot(data_cast_dna) +
  geom_point(aes(`TP53 +/-`, `DNASE1L3 -/-`), pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  stat_regline_equation(label.y = -0.5, aes(`TP53 +/-`, `DNASE1L3 -/-`, label = ..rr.label..)) +
  xlab("TP53 (+/-) vs HBC") + 
  ylab("DNASE1L3 (-/-) vs HBC") +
  ggtitle("Homozygous") + 
  theme
fig_reg1

data_cast_het <- reshape2::dcast(data_melt_het, motif ~ variable)
fig_reg2 <- ggplot(data_cast_het) +
  geom_point(aes(`TP53 +/-`, `DNASE1L3 +/-`), pch = 16, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  stat_regline_equation(label.y = -0.2, aes(`TP53 +/-`, `DNASE1L3 +/-`, label = ..rr.label..)) +
  xlab("TP53 (+/-) vs HBC") + 
  ylab("DNASE1L3 (+/-) vs HBC") +
  ggtitle("Heterozygous") + 
  theme
fig_reg2

figure <- ggarrange(fig_reg1, fig_reg2)
figure
ggsave(file.path(outdir, "fragment_end_contexts_dnase1l3_reg.pdf"), figure, width = 8, height = 4)
