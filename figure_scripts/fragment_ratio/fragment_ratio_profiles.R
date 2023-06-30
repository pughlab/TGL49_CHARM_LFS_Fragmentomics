library(tidyverse)
library(dplyr)
library(matrixStats)
library(grid)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_ratio"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_ratio"

### Import data (Starting with the 5Mb ratios)
data_ratio <- read.delim(file.path(path, "fragmentomics", "CHARM_LFS_ratio_5Mb.txt"))
data_normal <- read.delim(file.path(healthy, "fragmentomics", "TGL49_HBC_ratio_5Mb.txt"))
samples <- read.delim(file.path(path, "samples/sample_list.txt"))
de_table <- read.delim(file.path(outdir, "fragment_ratio_differential.txt"))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude | samples$germline_mutation == "unknown"), ]

### Keep only plasma samples and order based on clinical information
samples <- samples[samples$sWGS %in% colnames(data_ratio), ]
samples$diag <- "LFS"

### Format ratios
data_ratio <- data_ratio[, c("seqnames", "arm", "start", "end", samples$sWGS)]
data_ratio <- merge(data_ratio, data_normal, by = c("seqnames", "arm", "start", "end"))
data_ratio <- data_ratio[order(factor(data_ratio$seqnames, levels = c(paste0("chr", 1:22))),
                               data_ratio$start), ]
data_ratio$bin <- c(1:nrow(data_ratio))

data_melt <- reshape2::melt(data_ratio, id = c("seqnames", "arm", "start", "end", "bin"))
data_melt <- merge(data_melt, samples, by.x = "variable", by.y = "sWGS", all = TRUE)
data_melt[is.na(data_melt)] <- "HBC"
data_melt$diag <- ifelse(data_melt$cancer_status == "positive", "LFS-AC",
                         ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "yes", "LFS-PC",
                                ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "no", "LFS-H", "Healthy")))
data_melt$diag <- factor(data_melt$diag, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))
data_melt$mutation_type <- factor(data_melt$mutation_type, levels = c("HBC", "1", "2", "3", "5", "LOF", "Splice"),
                                  labels = c("Healthy", "Missense 1", "Missense 2", "Missense 3", "Missense 5", "LOF", "Splice"))
data_melt <- data_melt[order(data_melt$variable,
                             data_melt$bin), ]

### Make medians for each cluster
rows <- data_ratio[, c(1:4)]
rows$bin <- c(1:nrow(rows))
cohort <- rowMedians(as.matrix(data_ratio[, colnames(data_ratio) %in% samples$sWGS]))
c2 <- rowMedians(as.matrix(data_ratio[, colnames(data_ratio) %in% samples$sWGS[samples$mutation_type == "2"]]))
c3 <- rowMedians(as.matrix(data_ratio[, colnames(data_ratio) %in% samples$sWGS[samples$mutation_type == "3"]]))
c5 <- rowMedians(as.matrix(data_ratio[, colnames(data_ratio) %in% samples$sWGS[samples$mutation_type == "5"]]))
lof <- rowMedians(as.matrix(data_ratio[, colnames(data_ratio) %in% samples$sWGS[samples$mutation_type == "LOF"]]))
splice <- rowMedians(as.matrix(data_ratio[, colnames(data_ratio) %in% samples$sWGS[samples$mutation_type == "Splice"]]))

median_table <- cbind(rows, cohort, c2, c3, c5, lof, splice)
median_melt <- reshape2::melt(median_table, id = c("seqnames", "arm" , "start", "end", "bin"))
median_melt$variable <- factor(median_melt$variable, levels = c("cohort", "c2", "c3", "c5", "lof", "splice"),
                               labels = c("Cohort", "Missense 2", "Missense 3", "Missense 5", "LOF", "Splice"))
median_melt <- median_melt[order(median_melt$variable,
                             median_melt$bin), ]

### Make healthy median
median <- rowMedians(as.matrix(data_normal[, !(colnames(data_normal) %in% colnames(rows))]))
sd <- rowSds(as.matrix(data_normal[, !(colnames(data_normal) %in% colnames(rows))]))
hbc_median <- cbind(rows, median, sd)

### Compare variability between LFS and healthy
lfs_sd <- rowSds(as.matrix(data_ratio[, colnames(data_ratio)[colnames(data_ratio) %in% samples$sWGS]]))
neg_sd <- rowSds(as.matrix(data_ratio[, colnames(data_ratio)[colnames(data_ratio) %in% samples$sWGS[samples$cancer_status == "negative" & samples$previous_cancer == "no"]]]))
pc_sd <- rowSds(as.matrix(data_ratio[, colnames(data_ratio)[colnames(data_ratio) %in% samples$sWGS[samples$cancer_status == "negative" & samples$previous_cancer == "yes"]]]))
pos_sd <- rowSds(as.matrix(data_ratio[, colnames(data_ratio)[colnames(data_ratio) %in% samples$sWGS[samples$cancer_status == "positive"]]]))
diff_lfs <- mean(lfs_sd/hbc_median$sd)
diff_neg <- mean(neg_sd/hbc_median$sd)
diff_pc <- mean(pc_sd/hbc_median$sd)
diff_pos <- mean(pos_sd/hbc_median$sd)

# Plot profiles
mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line = element_line(colour = 'black', size = 0.5),
  strip.text.x = element_text(size=8),
  strip.text.y = element_text(size=8),
  axis.title.x = element_text(size=10),
  axis.title.y = element_text(size=10),
  axis.text.y = element_text(size=8),
  plot.title = element_text(size=14),
  legend.position = "none",
  legend.title = element_text(size=10),
  legend.text = element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill="white", color="white"),
  panel.spacing.x = unit(0.1, "lines"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
data_melt$arm <- factor(data_ratio$arm, levels=armlevels)
median_melt$arm <- factor(data_ratio$arm, levels=armlevels)
hbc_median$arm <- factor(data_ratio$arm, levels=armlevels)

arm <- data_ratio %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

# Generate Fragmentation plots
#Mutation Type
g1 <- ggplot() + 
  geom_line(data = data_melt, aes(x = bin, y = value, group = variable), linewidth = 0.5, alpha = 0.25) + 
  ggtitle("Fragmentation Profiles") +
  labs(x="Chromosome", y="Normalized Fragment Ratio", color="") + 
  facet_grid(mutation_type ~ arm, switch = "x",space = "free_x", scales = "free_x", labeller = labeller(arm = arm.labels)) + 
  coord_cartesian(xlim = NULL, ylim=c(-.08,.08), expand = TRUE) + 
  mytheme
g1

ggsave(file.path(outdir, "fragment_genome_profiles_type.pdf"), g1, width = 12, height = 5)

#LFS vs HBC
g2 <- ggplot() + 
  geom_line(data = data_melt, aes(x = bin, y = value, group = variable), linewidth = 0.25, alpha = 0.25) + 
  ggtitle("Fragmentation Profiles") +
  labs(x="Chromosome", y="Normalized Fragment Ratio", color="") + 
  facet_grid(diag ~ arm, switch = "x",space = "free_x", scales = "free_x", labeller = labeller(arm = arm.labels)) + 
  coord_cartesian(xlim = NULL, ylim=c(-.08,.08), expand = TRUE) + 
  mytheme
g2

ggsave(file.path(outdir, "fragment_genome_profiles_status.pdf"), g2, width = 12, height = 4)

#Medians
g3 <- ggplot() + 
  geom_line(data = hbc_median, aes(x = bin, y = median), color = "black", linewidth = 0.5, alpha = 0.75) +
  geom_ribbon(data = hbc_median, aes(x = bin, ymin = median - sd, ymax = median + sd), fill = "black", alpha = 0.25) +
  geom_line(data = median_melt, aes(x = bin, y = value), color = "red", linewidth = 0.5, alpha = 0.5) + 
  ggtitle("Fragmentation Profiles") +
  labs(x="Chromosome", y="Normalized Fragment Ratio", color="") + 
  guides(color = guide_legend(nrow = 1)) +
  facet_grid(variable ~ arm, switch = "x",space = "free_x", scales = "free_x", labeller = labeller(arm = arm.labels)) + 
  coord_cartesian(xlim = NULL, ylim=c(-.05,.05), expand = TRUE) + 
  mytheme
g3

ggsave(file.path(outdir, "fragment_genome_profiles_medians.pdf"), g3, width = 12, height = 4.5)

