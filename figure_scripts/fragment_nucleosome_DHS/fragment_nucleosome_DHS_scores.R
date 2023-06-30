library(tidyverse)
library(dplyr)
library(reshape2)
library(ggformula)
library(signal)
library(matrixStats)
library(lemon)
library(ggpubr)

### Set paths
path <- "data/griffin/DHS"
outdir <- ""
healthy_path <- "hbc/griffin/DHS"

### Read in samples
data_samples <- read.delim("sample_list.txt")
healthy_samples <- read.delim("hbc/sample_list.txt")

### Get list of files
data_griffin <- list.files(path, "features", full.names = TRUE)
data_normal <- list.files(healthy_path, "features", full.names = TRUE)

### Import data 
datalist <- lapply(data_griffin, function(x){read.delim(file = x)})
data_griffin <- do.call(rbind, datalist)
datalist <- lapply(data_normal, function(x){read.delim(file = x)})
data_normal <- do.call(rbind, datalist)

### Failed samples and remove
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG", "TGL49_0019_Cf_U_PE_317_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude, abnormal)),]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_griffin), ]

### Format griffin data
row.names(data_griffin) <- data_griffin$features
data_griffin <- data_griffin[, colnames(data_griffin) %in% data_samples$sWGS]
data_griffin <- as.data.frame(t(data_griffin))

row.names(data_normal) <- data_normal$features
data_normal <- data_normal[, -1]
data_normal <- as.data.frame(t(data_normal))

### Merge data and sample lists
data <- bind_rows(data_griffin, data_normal)
data_samples$diag <- "LFS"
samples <- as.data.frame(row.names(data_normal))
colnames(samples) <- "sWGS"
samples <- bind_rows(data_samples, samples)
samples[is.na(samples)] <- "HBC"

### Scale midpoint by mean coverage
data_mid <- data[, colnames(data) %like% "midpoint"]
data_cov <- data[, colnames(data) %like% "coverage"]
data_amp <- data[, colnames(data) %like% "amp"]
data_mid <- (data_mid - data_cov) + 1
data <- bind_cols(data_mid, data_amp)

### Melt data
data$sample <- row.names(data)

data_melt <- reshape2::melt(data, id = "sample")
data_melt$feature <- gsub("(.*_\\s*(.*$))", "\\2", data_melt$variable)
data_melt$site <- sub("_[^_]+$", "", data_melt$variable)
data_melt <- merge(data_melt, samples, by.x = "sample", by.y = "sWGS")

### Factor melt data
data_melt$cancer_status <- factor(data_melt$cancer_status, levels = c("HBC", "negative", "positive"),
                                  labels = c("Healthy", "Negative", "Positive"))

data_melt <- data_melt[data_melt$site %in% c("Cardiac", "Digestive", "Lymphoid", "Musculoskeletal", "Myeloid_erythroid", 
                                             "Neural", "Pulmonary", "Stromal_A", "Stromal_B", "Vascular_endothelial"), ]
data_melt$name <- factor(data_melt$site, levels = c("Lymphoid", "Myeloid_erythroid", "Vascular_endothelial", "Musculoskeletal", 
                                                    "Pulmonary", "Cardiac", "Digestive", "Neural", "Stromal_A", "Stromal_B"),
                         labels = c("Lymphoid", "Myeloid", "Endothelial", "Musculo\nSkeletal", "Pulmonary", "Cardiac", "Digestive", "Neural", "Stromal A", "Stromal B"))

### Calculate stats
data_stats <- data_melt %>%
  group_by(feature, site, name, cancer_status) %>%
  dplyr::summarise(mean = mean(value),
            median = median(value),
            sd = sd(value),
            N = n())
data_stats$variable <- paste0(data_stats$site, "_", data_stats$feature)

t_test <- c()
vars <- unique(data_stats$variable)
for (var in vars) {
  a <- t.test(data_melt$value[data_melt$cancer_status == "Healthy" & data_melt$variable == var], 
              data_melt$value[data_melt$cancer_status == "Negative" & data_melt$variable == var])$p.value
  b <- t.test(data_melt$value[data_melt$cancer_status == "Healthy" & data_melt$variable == var], 
              data_melt$value[data_melt$cancer_status == "Positive" & data_melt$variable == var])$p.value
  t_test <- c(t_test, 1, a, b)
}
t_test <- p.adjust(t_test, method = "BH")
data_stats$pvalue <- t_test
data_stats$annot <- ifelse(data_stats$pvalue < 0.05 & data_stats$pvalue > 0.01, "*",
                               ifelse(data_stats$pvalue < 0.01 & data_stats$pvalue > 0.001, "**",
                                      ifelse(data_stats$pvalue < 0.001, "***", "")))

### Format plotting tables
data_melt$feature <- factor(data_melt$feature, levels = c("midpoint", "coverage", "amp"),
                            labels = c("Midpoint\nCoverage", "Mean\nCoverage", "Amplitude"))
data_stats$feature <- factor(data_stats$feature, levels = c("midpoint", "coverage", "amp"),
                             labels = c("Midpoint\nCoverage", "Mean\nCoverage", "Amplitude"))
data_stats$place <- ifelse(data_stats$feature == "Midpoint\nCoverage", 1.15, 0.4)

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
               axis.title = element_text(size = 12),
               axis.line = element_line(colour = "black"),
               axis.text.x = element_text(size = 0, angle = 90, hjust = 1, vjust = 0.5),
               axis.text = element_text(size = 10),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               panel.spacing.x = unit(-1, "lines"),
               legend.position = "bottom",
               legend.background = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.key = element_rect(fill = "white"),
               legend.spacing.y = unit(0, "mm"),
               legend.key.size = unit(3, "mm"),
               strip.placement = "outside",
               strip.background = element_rect(color = NA, fill = NA),
               strip.text = element_text(face = "bold", size = 8))

### Graph Figure
data_melt <- data_melt %>%
  group_by(diag, feature, site) %>%
  arrange(value)
site_plot <- ggplot(data_melt, aes(cancer_status, value)) +
  geom_point(aes(group = cancer_status, color = cancer_status), alpha = 0.75,
             position = position_dodge2(.75), size = 0.5, pch = 16) +
  geom_text(data = data_stats, aes(x = cancer_status, y = place, label = annot), size = 4) +
  facet_grid(feature~name, scales = "free", switch = "y") +
  facet_rep_grid(feature~name, scales = "free", switch = "y") +
  labs(color = "Diagnosis") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  xlab("") + 
  ylab("") +
  scale_color_manual(values = c("grey65", "black", "red")) +
  theme +
  stat_summary(fun = match.fun(median), geom = "crossbar", size = 0.25, width = 0.75, group = 1, show.legend = FALSE, color = "#FF0000")
site_plot

ggsave(file.path(outdir, paste0("fragment_nucleosome_DHS_scores.pdf")), site_plot, device = "pdf", width = 8, height = 4, units = "in", limitsize = FALSE)

