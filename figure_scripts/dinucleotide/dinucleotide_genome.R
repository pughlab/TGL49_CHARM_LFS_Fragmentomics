library(gsignal)
library(dplyr)
library(ggpubr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(pracma)
library(ggh4x)
library(zplyr)

### Set paths
path <- "data/dinucleotide"
outdir <- ""
healthy_path <- "hbc/dinucleotide"
sim_data <- "data/DSP_sim_nucfreq.txt"

### Find paths
data <- list.files(path, "dinucleotide.txt", full.names = TRUE)
data <- data[grepl("genome", data)]
data_normal <- list.files(healthy_path, "dinucleotide.txt", full.names = TRUE)
data_normal <- data_normal[grepl("genome", data_normal)]

### Import data 
data <- read.delim(data)
data_normal <- read.delim(data_normal)
data_samples <- read.delim("sample_list.txt")
data_sim <- read.delim(sim_data)

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
#abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]
data <- data[data$sample %in% data_samples$sWGS, ]

### Sum frequencies across A/T and G/C combinations
data$context <- ifelse(data$context %in% c("AA", "AT", "TA", "TT"), "at",
                       ifelse(data$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
data <- data[complete.cases(data), ]
data_contexts <- aggregate(.~context+sample, data, sum)

data_normal$context <- ifelse(data_normal$context %in% c("AA", "AT", "TA", "TT"), "at",
                              ifelse(data_normal$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
data_normal <- data_normal[complete.cases(data_normal), ]
data_normal_contexts <- aggregate(.~context+sample, data_normal, sum)

### Format simulated data
data_sim <- data_sim[data_sim$Base %in% c("AA", "AT", "TA", "TT", "CC", "CG", "GC", "GG"), ]
data_sim$Base <- ifelse(data_sim$Base %in% c("AA", "AT", "TA", "TT"), "at",
                        ifelse(data_sim$Base %in% c("CC", "CG", "GC", "GG"), "gc", NA))
data_sim <- data_sim[complete.cases(data_sim), !(colnames(data_sim) %in% "Type")]
data_sim <- aggregate(.~Base, data_sim, sum)

sim_at <- data_sim[data_sim$Base == "at", !(colnames(data_sim) == "Base")]
sim_at <- as.numeric(sim_at)
sim_gc <- data_sim[data_sim$Base == "gc", !(colnames(data_sim) == "Base")]
sim_gc <- as.numeric(sim_gc)

### Calculate Log2ratios
columns <- paste0("X", 1:266)

data_at <- data_contexts[data_contexts$context == "at", ]
row.names(data_at) <- data_at$sample
data_at <- data_at[, colnames(data_at) %in% columns]
data_at <- sweep(data_at, 2, sim_at, "/")
data_at <- log2(data_at)

data_gc <- data_contexts[data_contexts$context == "gc", ]
row.names(data_gc) <- data_gc$sample
data_gc <- data_gc[, colnames(data_gc) %in% columns]
data_gc <- sweep(data_gc, 2, sim_gc, "/")
data_gc <- log2(data_gc)

normal_at <- data_normal_contexts[data_normal_contexts$context == "at", ]
row.names(normal_at) <- normal_at$sample
normal_at <- normal_at[, colnames(normal_at) %in% columns]
normal_at <- sweep(normal_at, 2, sim_at, "/")
normal_at <- log2(normal_at)

normal_gc <- data_normal_contexts[data_normal_contexts$context == "gc", ]
row.names(normal_gc) <- normal_gc$sample
normal_gc <- normal_gc[, colnames(normal_gc) %in% columns]
normal_gc <- sweep(normal_gc, 2, sim_gc, "/")
normal_gc <- log2(normal_gc)

### Format dataframes
distance <- c(-133:132)

data_at <- as.data.frame(t(data_at))
row.names(data_at) <- distance

data_gc <- as.data.frame(t(data_gc))
row.names(data_gc) <- distance

normal_at <- as.data.frame(t(normal_at))
row.names(normal_at) <- distance

normal_gc <- as.data.frame(t(normal_gc))
row.names(normal_gc) <- distance

### Calculate the healthy median
normal_at_median <- rowMedians(as.matrix(normal_at))
normal_at_sd <- rowSds(as.matrix(normal_at))

normal_gc_median <- rowMedians(as.matrix(normal_gc))
normal_gc_sd <- rowSds(as.matrix(normal_gc))

### Calculate the LFS median
at_median <- rowMedians(as.matrix(data_at))
at_sd <- rowSds(as.matrix(data_at))

gc_median <- rowMedians(as.matrix(data_gc))
gc_sd <- rowSds(as.matrix(data_gc))

samples_neg <- data_samples[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no", ]
data_at_neg <- data_at[, colnames(data_at) %in% samples_neg$sWGS]
at_median_neg <- rowMedians(as.matrix(data_at_neg))
at_sd_neg <- rowSds(as.matrix(data_at_neg))

data_gc_neg <- data_gc[, colnames(data_gc) %in% samples_neg$sWGS]
gc_median_neg <- rowMedians(as.matrix(data_gc_neg))
gc_sd_neg <- rowSds(as.matrix(data_gc_neg))

samples_sur <- data_samples[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes", ]
data_at_sur <- data_at[, colnames(data_at) %in% samples_sur$sWGS]
at_median_sur <- rowMedians(as.matrix(data_at_sur))
at_sd_sur <- rowSds(as.matrix(data_at_sur))

data_gc_sur <- data_gc[, colnames(data_gc) %in% samples_sur$sWGS]
gc_median_sur <- rowMedians(as.matrix(data_gc_sur))
gc_sd_sur <- rowSds(as.matrix(data_gc_sur))

samples_pos <- data_samples[data_samples$cancer_status == "positive", ]
data_at_pos <- data_at[, colnames(data_at) %in% samples_pos$sWGS]
at_median_pos <- rowMedians(as.matrix(data_at_pos))
at_sd_pos <- rowSds(as.matrix(data_at_pos))

data_gc_pos <- data_gc[, colnames(data_gc) %in% samples_pos$sWGS]
gc_median_pos <- rowMedians(as.matrix(data_gc_pos))
gc_sd_pos <- rowSds(as.matrix(data_gc_pos))

### Calculate z-scores
z_at <- (at_median - normal_at_median)/(sqrt(at_sd^2 + normal_at_sd^2))
z_gc <- (gc_median - normal_gc_median)/(sqrt(gc_sd^2 + normal_gc_sd^2))

z_at_neg <- (at_median_neg - normal_at_median)/(sqrt(at_sd_neg^2 + normal_at_sd^2))
z_gc_neg <- (gc_median_neg - normal_gc_median)/(sqrt(gc_sd_neg^2 + normal_gc_sd^2))

z_at_sur <- (at_median_sur - normal_at_median)/(sqrt(at_sd_sur^2 + normal_at_sd^2))
z_gc_sur <- (gc_median_sur - normal_gc_median)/(sqrt(gc_sd_sur^2 + normal_gc_sd^2))

z_at_pos <- (at_median_pos - normal_at_median)/(sqrt(at_sd_pos^2 + normal_at_sd^2))
z_gc_pos <- (gc_median_pos - normal_gc_median)/(sqrt(gc_sd_pos^2 + normal_gc_sd^2))

### Make median tables and format
normal_plot <- rbind(data.frame(distance = distance, median = normal_at_median, sd = normal_at_sd, context = "AA/AT/TA/TT", diag = "Healthy", test = "Log2(Observed/Expected)"),
                     data.frame(distance = distance, median = normal_gc_median, sd = normal_gc_sd, context = "CC/CG/GC/GG", diag = "Healthy", test = "Log2(Observed/Expected)"))
normal_plot$distance <- as.numeric(normal_plot$distance)

lfs_plot <- rbind(data.frame(distance = distance, median = at_median_neg, sd = at_sd_neg, context = "AA/AT/TA/TT", diag = "LFS-H", test = "Log2(Observed/Expected)"),
                  data.frame(distance = distance, median = gc_median_neg, sd = gc_sd_neg, context = "CC/CG/GC/GG", diag = "LFS-H", test = "Log2(Observed/Expected)"),
                  data.frame(distance = distance, median = at_median_sur, sd = at_sd_sur, context = "AA/AT/TA/TT", diag = "LFS-PC", test = "Log2(Observed/Expected)"),
                  data.frame(distance = distance, median = gc_median_sur, sd = gc_sd_sur, context = "CC/CG/GC/GG", diag = "LFS-PC", test = "Log2(Observed/Expected)"),
                  data.frame(distance = distance, median = at_median_pos, sd = at_sd_pos, context = "AA/AT/TA/TT", diag = "LFS-AC", test = "Log2(Observed/Expected)"),
                  data.frame(distance = distance, median = gc_median_pos, sd = gc_sd_pos, context = "CC/CG/GC/GG", diag = "LFS-AC", test = "Log2(Observed/Expected)"))
lfs_plot$distance <- as.numeric(lfs_plot$distance)

z_plot <- rbind(data.frame(distance = distance, median = z_at_neg, sd = 0, context = "AA/AT/TA/TT", diag = "LFS-H", test = "Z-Score"),
                data.frame(distance = distance, median = z_gc_neg, sd = 0, context = "CC/CG/GC/GG", diag = "LFS-H", test = "Z-Score"),
                data.frame(distance = distance, median = z_at_sur, sd = 0, context = "AA/AT/TA/TT", diag = "LFS-PC", test = "Z-Score"),
                data.frame(distance = distance, median = z_gc_sur, sd = 0, context = "CC/CG/GC/GG", diag = "LFS-PC", test = "Z-Score"),
                data.frame(distance = distance, median = z_at_pos, sd = 0, context = "AA/AT/TA/TT", diag = "LFS-AC", test = "Z-Score"),
                data.frame(distance = distance, median = z_gc_pos, sd = 0, context = "CC/CG/GC/GG", diag = "LFS-AC", test = "Z-Score"))
z_plot$distance <- as.numeric(z_plot$distance)

data_plot <- bind_rows(normal_plot, lfs_plot, z_plot)
data_plot$median <- as.numeric(data_plot$median)
data_plot$sd <- as.numeric(data_plot$sd)

data_plot$diag <- factor(data_plot$diag, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))

### Make regression tables and format
reg_plot <- rbind(data.frame(var1 = normal_at_median, var2 = at_median_neg, context = "AA/AT/TA/TT", diag = "LFS-H"),
                  data.frame(var1 = normal_gc_median, var2 = gc_median_neg, context = "CC/CG/GC/GG", diag = "LFS-H"),
                  data.frame(var1 = normal_at_median, var2 = at_median_sur, context = "AA/AT/TA/TT", diag = "LFS-PC"),
                  data.frame(var1 = normal_gc_median, var2 = gc_median_sur, context = "CC/CG/GC/GG", diag = "LFS-PC"),
                  data.frame(var1 = normal_at_median, var2 = at_median_pos, context = "AA/AT/TA/TT", diag = "LFS-AC"),
                  data.frame(var1 = normal_gc_median, var2 = gc_median_pos, context = "CC/CG/GC/GG", diag = "LFS-AC"))
reg_plot$var1 <- as.numeric(reg_plot$var1)
reg_plot$var2 <- as.numeric(reg_plot$var2)

reg_plot$diag <- factor(reg_plot$diag, levels = c("LFS-H", "LFS-PC", "LFS-AC"))

### Calculate Z-scores
data_zscore <- data_plot[data_plot$test == "Z-Score", ]
data_zscore$position <- ifelse(data_zscore$distance > 83 |
                                 data_zscore$distance < -83, "Flanking", "Nucleosome")
contexts <- unique(data_zscore$context) 
diags <- unique(data_zscore$diag)
positions <- unique(data_zscore$position)
x <- c()

for (context in contexts) {
  for (diag in diags) {
    a <- t.test(data_zscore$median[data_zscore$context == context & data_zscore$diag == diag & data_zscore$position == "Flanking"],
                data_zscore$median[data_zscore$context == context & data_zscore$diag == diag & data_zscore$position == "Nucleosome"])$p.value
    a <- formatC(a, format = "e", digits = 2)
    x <- c(x,a)
  }
}
a <- ks.test(data_zscore$median[data_zscore$context == "AA/AT/TA/TT" & data_zscore$position == "Flanking"],
             data_zscore$median[data_zscore$context == "AA/AT/TA/TT" & data_zscore$position == "Nucleosome"])$p.value
a <- formatC(a, format = "e", digits = 2)
b <- ks.test(data_zscore$median[data_zscore$context == "CC/CG/GC/GG" & data_zscore$position == "Flanking"],
             data_zscore$median[data_zscore$context == "CC/CG/GC/GG" & data_zscore$position == "Nucleosome"])$p.value
b <- formatC(b, format = "e", digits = 2)


ann_text <- data.frame(x = c(rep(0.2, 3), rep(0.55, 3)),
                       y = c(rep(0.55, 3), rep(0.55, 3)),
                       context = c(rep("AA/AT/TA/TT", 3), rep("CC/CG/GC/GG", 3)),
                       diag = rep(diags, 2),
                       pvalue = x,
                       label = c(paste0("KS-test\np-value\n", x)))

data_zscore$context <- factor(data_zscore$context, levels = c("CC/CG/GC/GG", "AA/AT/TA/TT"))
ann_text$context <- factor(ann_text$context, levels = c("CC/CG/GC/GG", "AA/AT/TA/TT"))

### Set plotting limits
min <- min(data_plot$median[data_plot$diag %in% c("Healthy", "LFS")]) - 0.025
max <- max(data_plot$median[data_plot$diag %in% c("Healthy", "LFS")]) + 0.025

min2 <- min(data_plot$median[data_plot$test == "Z-Score"]) - 0.025
max2 <- max(data_plot$median[data_plot$test == "Z-Score"]) + 0.025

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
               axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

### Plot curves
normal_plot <- normal_plot[,-5]
fig <- ggplot(data_plot[!(data_plot$diag == "Healthy"), ]) +
  geom_line(data = normal_plot, aes(distance, median, group = context, color = "aaaa"), size = 0.75) +
  geom_line(aes(distance, median, group = context, color = context), size = 0.75) +
  #geom_ribbon(aes(distance, ymin = median - sd, ymax = median + sd, group = context, fill = context), alpha = 0.5) +
  geom_vline(xintercept = c(-83, 83), linetype = "dashed", size = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  xlab("Position Relative to 167bp") + 
  ylab("") +
  labs(color = "Dinucleotide", fill = "Dinucleotide") +
  scale_color_manual(labels = c("Healthy", "AA/AT/TA/TT", "CC/CG/GC/GG"), 
                     values = c("aaaa" = "black", "AA/AT/TA/TT" = "#FF9333", "CC/CG/GC/GG" = "#339FFF")) +
  #scale_fill_manual(values = c("#FF9333", "#339FFF")) +
  ggtitle(paste0("Dinucleotide Frequencies")) + 
  facet_grid(diag~test, scales = "free_y") +
  facetted_pos_scales(y = list(scale_y_continuous(limits = c(min, max)),
                               scale_y_continuous(limits = c(min, max)),
                               scale_y_continuous(limits = c(min, max)),
                               scale_y_continuous(limits = c(min2, max2)),
                               scale_y_continuous(limits = c(min2, max2)),
                               scale_y_continuous(limits = c(min2, max2)))) +
  theme +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-132, 132), expand = c(0,0))
fig

ggsave(file.path(outdir, "dinucleotide_contexts_genome.pdf"), fig, width = 7.5, height = 6.5)

fig_reg <- ggplot(reg_plot, aes(var1, var2)) +
  geom_point(pch = 16, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  xlab("Healthy Control") + 
  ylab("") +
  ggtitle(paste0("Correlation to Healthy")) + 
  facet_grid(diag~context) +
  theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
fig_reg

ggsave(file.path(outdir, "dinucleotide_contexts_regression.pdf"), fig_reg, width = 5, height = 4)

fig_z <- ggplot(data_zscore) +
  geom_density(aes(median, fill = position), color = NA, alpha = 0.25) +
  geom_abs_text(data = ann_text, aes(xpos = x, ypos = y, label = label), size = 3.5) +
  xlab("Median Z-scores") + 
  ylab("Density") +
  labs(fill = "") +
  scale_fill_manual(values = c("red", "blue")) +
  ggtitle("LFS vs Healthy") + 
  facet_grid(diag ~ context, scales = "free", space = "fixed") +
  theme +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = c(-0.2, 0, 0.2))
fig_z
ggsave(file.path(outdir, paste0("dinucleotide_contexts_zscores.pdf")), fig_z, width = 6, height = 4)



