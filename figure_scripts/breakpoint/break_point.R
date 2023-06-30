library(tidyverse)
library(dplyr)
library(matrixStats)
library(reshape2)

### Set paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/breakpoint"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/breakpoint"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/breakpoint"

### Find paths
data <- list.files(path, "ratio.txt", full.names = TRUE)
data <- data[data %like% "CHARM_LFS"]
data_normal <- list.files(healthy_path, "ratio.txt", full.names = TRUE)

### Import data 
data <- read.delim(data)
data_normal <- read.delim(data_normal)
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
#abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]
data_samples <- data_samples[data_samples$sWGS %in% data$sample, ]

### Format healthy controls and find median
colnames(data_normal) <- c(-15:14, "nucleotide", "sample")
melt_normal <- reshape2::melt(data_normal, id = c("nucleotide", "sample"))

median_normal <- melt_normal %>%
  group_by(nucleotide, variable) %>%
  dplyr::summarise(median = median(value),
                   sd = sd(value))
median_normal <- median_normal[median_normal$variable %in% c(-14:14), ]
median_normal$status <- "Healthy"

### Format LFS data and find medians
samples_h <- data_samples[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no", ]
samples_p <- data_samples[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes", ]
samples_a <- data_samples[data_samples$cancer_status == "positive", ]

colnames(data) <- c(-15:14, "nucleotide", "sample")
melt_data <- reshape2::melt(data, id = c("nucleotide", "sample"))

melt_h <- melt_data[melt_data$sample %in% samples_h$sWGS, ]
median_h <- melt_h %>%
  group_by(nucleotide, variable) %>%
  dplyr::summarise(median = median(value),
                   sd = sd(value))
median_h <- median_h[median_h$variable %in% c(-14:14), ]
median_h$status <- "LFS-H"

melt_p <- melt_data[melt_data$sample %in% samples_p$sWGS, ]
median_p <- melt_p %>%
  group_by(nucleotide, variable) %>%
  dplyr::summarise(median = median(value),
                   sd = sd(value))
median_p <- median_p[median_p$variable %in% c(-14:14), ]
median_p$status <- "LFS-PC"

melt_a <- melt_data[melt_data$sample %in% samples_a$sWGS, ]
median_a <- melt_a %>%
  group_by(nucleotide, variable) %>%
  dplyr::summarise(median = median(value),
                   sd = sd(value))
median_a <- median_a[median_a$variable %in% c(-14:14), ]
median_a$status <- "LFS-AC"

### Combine plotting dataframe
data_plot <- bind_rows(median_normal, median_h, median_p, median_a)
data_plot$nucleotide <- factor(data_plot$nucleotide, levels = c("A", "T", "C", "G"))
data_plot$status <- factor(data_plot$status, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))

### Calculate statistics
sd_difference <- median(data_plot$sd[data_plot$status == "Healthy" & data_plot$variable %in% c(-10:10)])
sd_difference_h <- median(data_plot$sd[data_plot$status == "LFS-H" & data_plot$variable %in% c(-10:10)])/sd_difference
sd_difference_p <- median(data_plot$sd[data_plot$status == "LFS-PC" & data_plot$variable %in% c(-10:10)])/sd_difference
sd_difference_a <- median(data_plot$sd[data_plot$status == "LFS-AC" & data_plot$variable %in% c(-10:10)])/sd_difference

data_stats <- data.frame(nucleotide = c(),
                         variable = c(),
                         status = c())

for (i in c("A", "C", "G", "T")) {
  for (j in c(-14:14)) {
    w <- melt_normal[melt_normal$nucleotide == i & melt_normal$variable == j, ]
    x <- melt_h[melt_h$nucleotide == i & melt_h$variable == j, ]
    y <- melt_p[melt_p$nucleotide == i & melt_p$variable == j, ]
    z <- melt_a[melt_a$nucleotide == i & melt_a$variable == j, ]
    
    a <- t.test(w$value, x$value)$p.value
    b <- t.test(w$value, y$value)$p.value
    c <- t.test(w$value, z$value)$p.value
    
    p <- data.frame(nucleotide = i,
                    variable = j,
                    status = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"),
                    pvalue = c(1, a, b, c))
    data_stats <- bind_rows(data_stats, p)
  }
}

data_plot <- merge(data_plot, data_stats, by = c("nucleotide", "variable", "status"))
data_plot$annot <- ifelse(data_plot$pvalue < 0.05 & data_plot$pvalue > 0.01, "<0.05",
                                ifelse(data_plot$pvalue < 0.01 & data_plot$pvalue > 0.001, "<0.01",
                                       ifelse(data_plot$pvalue < 0.001, "<0.001", "NS")))

### Make a dataframe with the difference between healthy and LFS
data_difference <- reshape2::dcast(data_plot, nucleotide + variable ~ status, value.var = "median")
data_difference$`LFS-H` <- data_difference$`LFS-H`/data_difference$Healthy
data_difference$`LFS-PC` <- data_difference$`LFS-PC`/data_difference$Healthy
data_difference$`LFS-AC` <- data_difference$`LFS-AC`/data_difference$Healthy
data_difference <- reshape2::melt(data_difference, id = c("nucleotide", "variable"))
colnames(data_difference) <- c("nucleotide", "position", "variable", "value")
data_difference <- data_difference[data_difference$variable %like% "LFS", ]

data_difference <- merge(data_difference, data_stats, by.x = c("nucleotide", "position", "variable"), by.y = c("nucleotide", "variable", "status"))
data_difference$annot <- ifelse(data_difference$pvalue < 0.05 & data_difference$pvalue > 0.01, "<0.05",
                          ifelse(data_difference$pvalue < 0.01 & data_difference$pvalue > 0.001, "<0.01",
                                 ifelse(data_difference$pvalue < 0.001, "<0.001", "NS")))

### Caclulate difference stats
data_difference$di <- ifelse(data_difference$nucleotide %in% c("A", "T"), "AT", "GC")
data_difference$change <- ifelse(data_difference$value > 1, "increase", "decrease")

model <- chisq.test(table(data_difference$change, data_difference$di))
stats_chisq_change <- model$p.value

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               legend.position = "top",
               legend.key = element_blank(),
               legend.background = element_blank(),
               axis.text = element_text(size = 11),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               strip.background = element_blank(),
               strip.text = element_text(size = 13))

### Plot comparisons
fig_context <- ggplot(data_plot) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 15, linetype = "dashed", size = 0.5) +
  geom_point(aes(variable, median, color = nucleotide, size = annot)) +
  geom_line(aes(variable, median, group = nucleotide, color = nucleotide)) +
  geom_errorbar(aes(variable, ymin = median - sd, ymax = median + sd, color = nucleotide), width = 0.5) +
  xlab("Position from 5' Fragment End") + 
  ylab("Frequency (Observed/Expected)") +
  labs(size = "p-value") +
  facet_grid(status~nucleotide) +
  scale_color_manual(values =c("G" = "#E69F00", "A" = "#009E73", "T" = "#D55E00", "C" = "#0072B2"),
                     guide = "none") +
  scale_size_manual(values = c("NS" = NA, "<0.05" = 0.5, "<0.01" = 1, "<0.001" = 1.5), guide = "legend") +
  ggtitle("Breakpoint Nucleotide Contexts") + 
  theme +
  scale_x_discrete(breaks = c(-10, -5, 0, 5, 10), labels = c(-10, -5, 0, 5, 10))
fig_context

ggsave(file.path(outdir, paste0("breakpoint_contexts.pdf")), fig_context, width = 8, height = 5)

fig_diff <- ggplot(data_difference) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 15, linetype = "dashed", size = 0.5) +
  geom_point(aes(position, value, color = nucleotide, size = annot)) +
  geom_line(aes(position, value, group = nucleotide, color = nucleotide)) +
  xlab("Position from 5' Fragment End") + 
  ylab("Median Difference") +
  labs(size = "p-value") +
  facet_grid(variable~nucleotide) +
  scale_color_manual(values = c("G" = "#E69F00", "A" = "#009E73", "T" = "#D55E00", "C" = "#0072B2"),
                     guide = "none") +
  ggtitle("Median Difference vs Healthy Controls") + 
  theme +
  scale_size_manual(values = c("NS" = NA, "<0.05" = 0.5, "<0.01" = 1, "<0.001" = 1.5), guide = "legend") +
  scale_x_discrete(breaks = c(-10, -5, 0, 5, 10), labels = c(-10, -5, 0, 5, 10))
fig_diff

ggsave(file.path(outdir, paste0("breakpoint_difference.pdf")), fig_diff, width = 8, height = 5)

