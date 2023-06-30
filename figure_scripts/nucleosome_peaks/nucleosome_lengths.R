library(dplyr)
library(matrixStats)
library(cowplot)
library(ggpubr)
library(gridExtra)

### Set paths
path <- "data/nucleosome_peaks"
outdir <- ""
healthy_path <- "hbc/nucleosome_peaks"
samples <- "sample_list.txt"

## Import data
data_prox <- read.delim(list.files(path, "proximal", full.names = TRUE))
data_dist <- read.delim(list.files(path, "distal", full.names = TRUE))
normal_prox <- read.delim(list.files(healthy_path, "proximal", full.names = TRUE))
normal_dist <- read.delim(list.files(healthy_path, "distal", full.names = TRUE))
data_samples <- read.delim(samples)

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[(data_samples$cancer_status %in% c("positive", "negative")), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_prox), ]

### Format data frames
data_prox <- data_prox[ , colnames(data_prox) %in% data_samples$sWGS]
names_prox <- colnames(data_prox)
data_dist <- data_dist[ , colnames(data_dist) %in% data_samples$sWGS]
names_dist <- colnames(data_dist)
normal_prox <- normal_prox[ , -1]
normal_dist <- normal_dist[ , -1]

### Calculate medians
data_prox_median <- rowMedians(as.matrix(data_prox))
data_dist_median <- rowMedians(as.matrix(data_dist))
normal_prox_median <- rowMedians(as.matrix(normal_prox))
normal_dist_median <- rowMedians(as.matrix(normal_dist))

### Make previvor median
samples_neg <- data_samples[data_samples$cancer_status == "negative" &
                            data_samples$previous_cancer == "no", ]
data_neg <- data_prox[ , colnames(data_prox) %in% samples_neg$sWGS]
data_prox_neg <- rowMedians(as.matrix(data_neg))
data_neg <- data_dist[ , colnames(data_dist) %in% samples_neg$sWGS]
data_dist_neg <- rowMedians(as.matrix(data_neg))

### Make cancer survivor median
samples_sur <- data_samples[data_samples$cancer_status == "negative" &
                              data_samples$previous_cancer == "yes", ]
data_sur <- data_prox[ , colnames(data_prox) %in% samples_sur$sWGS]
data_prox_sur <- rowMedians(as.matrix(data_sur))
data_sur <- data_dist[ , colnames(data_dist) %in% samples_sur$sWGS]
data_dist_sur <- rowMedians(as.matrix(data_sur))

### Make cancer positive median
samples_pos <- data_samples[data_samples$cancer_status == "positive", ]
data_pos <- data_prox[ , colnames(data_prox) %in% samples_pos$sWGS]
data_prox_pos <- rowMedians(as.matrix(data_pos))
data_pos <- data_dist[ , colnames(data_dist) %in% samples_pos$sWGS]
data_dist_pos <- rowMedians(as.matrix(data_pos))

### Make median table and format for plotting
data_plot1 <- rbind(data.frame(length = c(1:600),
                               proximal = data_prox_neg,
                               distal = data_dist_neg,
                               diag = "LFS-H"),
                    data.frame(length = c(1:600),
                               proximal = data_prox_sur,
                               distal = data_dist_sur,
                               diag = "LFS-PC"),
                    data.frame(length = c(1:600),
                               proximal = data_prox_pos,
                               distal = data_dist_pos,
                               diag = "LFS-AC"))
data_plot1 <- reshape2::melt(data_plot1, id = c("length", "diag"))
data_plot1$diag <- factor(data_plot1$diag, levels = c("LFS-H", "LFS-PC", "LFS-AC"))
data_plot1$variable <- factor(data_plot1$variable, levels = c("proximal", "distal"),
                              labels = c("Proximal (+/- 1kb)", "Distal (>1kb)"))

data_median <- data.frame(length = c(1:600),
                          proximal = normal_prox_median,
                          distal = normal_dist_median)                
data_median <- reshape2::melt(data_median, id = "length")
data_median$variable <- factor(data_median$variable, levels = c("proximal", "distal"),
                               labels = c("Proximal (+/- 1kb)", "Distal (>1kb)"))

### Calculate Proportions
prop_prox <- colSums(data_prox[row.names(data_prox) %in% c(10:150), ])/100
prop_dist <- colSums(data_dist[row.names(data_dist) %in% c(10:150), ])/100
normal_prop_prox <- colSums(normal_prox[row.names(normal_prox) %in% c(10:150), ])/100
normal_prop_dist <- colSums(normal_dist[row.names(normal_dist) %in% c(10:150), ])/100

### Make proportion table
data_plot2  <- merge(data.frame(sample = names_prox, proximal = prop_prox),
                     data.frame(sample = names_dist, distal = prop_dist), by = "sample")
data_plot2 <- merge(data_plot2, data_samples[,colnames(data_samples) %in% c("cancer_status", "previous_cancer", "sWGS")], by.x = "sample", by.y = "sWGS")
data_plot2 <- bind_rows(data_plot2,
                        data.frame(proximal = normal_prop_prox,
                                   distal = normal_prop_dist))
data_plot2[is.na(data_plot2)] <- "Healthy"
data_plot2$diag <- ifelse(data_plot2$cancer_status == "negative" & data_plot2$previous_cancer == "no", "LFS-H",
                          ifelse(data_plot2$cancer_status == "negative" & data_plot2$previous_cancer == "yes", "LFS-PC",
                                 ifelse(data_plot2$cancer_status == "positive", "LFS-AC", 
                                        ifelse(data_plot2$cancer_status == "Healthy", "Healthy", NA))))
data_plot2 <- reshape2::melt(data_plot2[,c(1:3,6)], id = c("sample", "diag"))

data_plot2$diag <- factor(data_plot2$diag, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))
data_plot2$variable <- factor(data_plot2$variable, levels = c("proximal", "distal"),
                              labels = c("Proximal (+/- 1kb)", "Distal (>1kb)"))

### Calculate stats
data_stats <- data_plot2 %>%
  group_by(variable, diag) %>%
  dplyr::summarise(mean=mean(value),
                   median=median(value),
                   sd=sd(value),
                   N=n())

variables <- unique(data_stats$variable)
diags <- unique(data_stats$diag)
x <- c()
for (variable in variables) {
  for (diag in diags) {
    a <- t.test(data_plot2$value[data_plot2$variable == variable & data_plot2$diag == "Healthy"],
                data_plot2$value[data_plot2$variable == variable & data_plot2$diag == diag])$p.value
    x <- c(x,a)
  }
}
data_stats$pvalue <- x
data_stats$annot <- ifelse(data_stats$pvalue < 0.05 & data_stats$pvalue > 0.01, "*",
                           ifelse(data_stats$pvalue < 0.01 & data_stats$pvalue > 0.001, "**",
                                  ifelse(data_stats$pvalue < 0.001, "***", "")))

data_aov <- data_plot2[data_plot2$variable == "Proximal (+/- 1kb)", ]
model <- aov(value ~ diag, data = data_aov)
x <- summary(model)
y <- round(x[[1]][["Pr(>F)"]], 4)
stats_anova_proximal <- y[[1]]

data_aov <- data_plot2[data_plot2$variable == "Distal (>1kb)", ]
model <- aov(value ~ diag, data = data_aov)
x <- summary(model)
y <- round(x[[1]][["Pr(>F)"]], 4)
stats_anova_distal <- y[[1]]

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
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(size = 13))

plot1 <- ggplot(data_plot1) +
  geom_line(aes(length, value, color = diag)) +
  geom_line(data = data_median, aes(length, value)) +
  geom_vline(xintercept = c(168), linetype = "dashed", color = "black") +
  facet_grid(variable~diag) +
  ggtitle("Fragment Lengths vs Distance from Nucleosome") +
  xlab("Fragment Length") + 
  ylab("Frequency (%)") +
  scale_color_manual(values =c("#1F78B4", "#33A02C", "#E31A1C")) +
  theme +
  scale_x_continuous(limits = c(90,220)) +
  scale_y_continuous(limits=c(0, 3), expand = c(0,0))
plot1

ggsave(file.path(outdir, "nucleosome_peaks_frequency.pdf"), width = 9, height = 3.5)

plot2 <- ggplot(data_plot2) +
  geom_boxplot(aes(diag, value, fill = diag), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(diag, value), width = 0.2, size = 0.75) +
  geom_text(data = data_stats, aes(diag, 0.39, label = annot), size = 4) +
  geom_text(data = data_stats, aes(diag, 0.06, label = N), size = 4) +
  facet_grid(.~variable) +
  ggtitle("Fragment Length Proportions") +
  labs(fill = "Patient Type") +
  xlab("Proportion (>150bp)") + 
  ylab("Patient Type") +
  scale_fill_manual(values =c("grey65", "#1F78B4", "#33A02C", "#E31A1C")) +
  theme +
  theme(axis.text.x = element_blank(),
        legend.position = "right") +
  scale_y_continuous(limits = c(0.05, 0.4))
plot2

ggsave(file.path(outdir, "nucleosome_peaks_proportions.pdf"), width = 6, height = 3.5)



