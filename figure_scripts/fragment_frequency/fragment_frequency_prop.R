library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(ggh4x)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/insert_size"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_frequency"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/insert_size"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"

## Import data
data_frequency <- read.delim(list.files(path, "LFS_fragment.txt", full.names = TRUE))
data_normal <- read.delim(list.files(healthy_path, "fragment.txt", full.names = TRUE))
data_samples <- read.delim(samples)

### Restrict size (10-500bp)
data_frequency <- data_frequency[data_frequency$length %in% c(10:500), ]
data_normal <- data_normal[data_normal$length %in% c(10:500), ]

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[(data_samples$cancer_status %in% c("positive", "negative")), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_frequency), ]
row.names(data_frequency) <- data_frequency$length
data_frequency <- data_frequency[ , colnames(data_frequency) %in% data_samples$sWGS]
data_frequency[is.na(data_frequency)] <- 0

row.names(data_normal) <- data_normal$length
data_normal <- data_normal[, !(colnames(data_normal) == "length")]

### Calculate proportions of size bins
frag_1 <- colSums(as.matrix(data_frequency[row.names(data_frequency) %in% c(10:150), ]))/colSums(data_frequency[row.names(data_frequency) %in% c(10:250), ])
frag_2 <- colSums(as.matrix(data_frequency[row.names(data_frequency) %in% c(151:180), ]))/colSums(data_frequency)
frag_3 <- colSums(as.matrix(data_frequency[row.names(data_frequency) %in% c(250:500), ]))/colSums(data_frequency)

norm_1 <- colSums(as.matrix(data_normal[row.names(data_normal) %in% c(10:150), ]))/colSums(data_normal[row.names(data_normal) %in% c(10:250), ])
norm_2 <- colSums(as.matrix(data_normal[row.names(data_normal) %in% c(151:180), ]))/colSums(data_normal)
norm_3 <- colSums(as.matrix(data_normal[row.names(data_normal) %in% c(250:500), ]))/colSums(data_normal)

### Make proportions table
data_frag <- data.frame(sample = colnames(data_frequency),
                        frag_1 = frag_1,
                        frag_2 = frag_2,
                        frag_3 = frag_3)
data_frag <- merge(data_frag, data_samples, by.x = "sample", by.y = "sWGS")
data_norm <- data.frame(sample = colnames(data_normal),
                        frag_1 = norm_1,
                        frag_2 = norm_2,
                        frag_3 = norm_3)
data <- bind_rows(data_frag, data_norm)
data[is.na(data)] <- "Healthy"
data$diag <- ifelse(data$cancer_status == "positive", "LFS-AC",
                    ifelse(data$cancer_status == "negative" & data$previous_cancer == "yes", "LFS-PC",
                           ifelse(data$cancer_status == "negative" & data$previous_cancer == "no", "LFS-H", "Healthy")))
data <- data[colnames(data) %in% c("sample", "frag_1", "frag_2", "frag_3", "diag")]
data_melt <- reshape2::melt(data, id = c("sample", "diag"))
data_melt$diag <- factor(data_melt$diag, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))
data_melt$variable <- factor(data_melt$variable, levels = c("frag_1", "frag_2", "frag_3"),
                             labels = c("10-150bp", "151-180bp", "250-500bp"))

### Calculate stats
data_stats <- data_melt %>%
  group_by(diag, variable) %>%
  dplyr::summarise(mean=mean(value),
                   sd=sd(value),
                   N=n())

patients <- unique(data_stats$diag)
sizes <- unique(data_stats$variable)
t_test <- c()
for (patient in patients) {
  for (size in sizes) {
    a <- t.test(data_melt$value[data_melt$diag == "Healthy" & data_melt$variable == size], data_melt$value[data_melt$diag == patient & data_melt$variable == size])$p.value
    t_test <- c(t_test, a)
  }
}
data_stats$pvalue <- t_test
data_stats$annot <- ifelse(data_stats$pvalue < 0.05 & data_stats$pvalue > 0.01, "*",
                           ifelse(data_stats$pvalue < 0.01 & data_stats$pvalue > 0.001, "**",
                                  ifelse(data_stats$pvalue < 0.001, "***", "")))

### Set limits
data_stats$max <- ifelse(data_stats$variable == "10-150bp", 0.41,
                         ifelse(data_stats$variable == "151-180bp", 0.63, 0.23))
data_stats$min <- ifelse(data_stats$variable == "10-150bp", 0.05,
                         ifelse(data_stats$variable == "151-180bp", 0.25, 0))


### Plot differences
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.title = element_text(size = 13),
               axis.line = element_line(colour = "black"),
               axis.text = element_text(size = 13),
               axis.text.x = element_blank(),
               legend.position = "none",
               legend.text = element_text(size = 13),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(size = 13))

plot <- ggplot(data_melt) +
  geom_boxplot(aes(diag, value, fill = diag), alpha = 0.5) +
  geom_text(data = data_stats, aes(diag, max, label = annot), size = 6) +
  geom_text(data = data_stats, aes(diag, min, label = N), size = 4) +
  ggtitle("Fragment Size Proportions") +
  xlab("Cancer Status") + 
  ylab("Proportion of Fragments") +
  facet_wrap(.~variable, scales = "free_y") +
  facetted_pos_scales(y = list(scale_y_continuous(limits = c(0.05, 0.43)),
                               scale_y_continuous(limits = c(0.25, 0.65)),
                               scale_y_continuous(limits = c(0, 0.25)))) +
  theme +
  scale_fill_manual(name = " ", values =c("black", "#1F78B4", "#33A02C", "#E31A1C"))
  
plot

ggsave(file.path(outdir, "fragment_proportions.pdf"), width = 4.5, height = 3)

