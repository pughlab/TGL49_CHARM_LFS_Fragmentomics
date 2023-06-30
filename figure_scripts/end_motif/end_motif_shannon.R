library(tidyverse)
library(dplyr)
library(matrixStats)
library(cowplot)
library(ggh4x)
library(zplyr)
library(patchwork)

### Set paths
path <- "data/end_motifs"
outdir <- ""
healthy_path <- "hbc/end_motifs"

### Find paths
data <- list.files(path, "motifs.txt", full.names = TRUE)
data <- data[grepl("genome2", data)]
data_normal <- list.files(healthy_path, "motifs.txt", full.names = TRUE)
data_normal <- data_normal[grepl("genome2", data_normal)]

### Import data 
data <- read.delim(data)
data_normal <- read.delim(data_normal)
data_samples <- read.delim("sample_list.txt")

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
#abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]

### Calculate Shannon Index
data_shan <- data[,-1]
data_shan <- colSums((-data_shan*log(data_shan))/log(nrow(data_shan)))
hbc_shan <- data_normal[,-1]
hbc_shan <- colSums((-hbc_shan*log(hbc_shan))/log(nrow(hbc_shan)))

### Calculate Gini Index
data_gini <- data[,-1]
data_gini <- 1 - colSums(data_gini^2)
hbc_gini <- data_normal[,-1]
hbc_gini <- 1 - colSums(hbc_gini^2)

### Combine data
data_score <- data.frame(shannon = data_shan,
                         gini = data_gini)
data_score$diag <- "LFS"
data_score$sample <- row.names(data_score)
data_score <- merge(data_score, data_samples, by.x = "sample", by.y = "sWGS")

hbc_score <- data.frame(shannon = hbc_shan,
                        gini = hbc_gini)
hbc_score$diag <- "HBC"
hbc_score$sample <- row.names(hbc_score)

data_plot <- bind_rows(data_score, hbc_score)
data_plot[is.na(data_plot)] <- "HBC"

data_plot$status <- ifelse(data_plot$cancer_status == "negative" & data_plot$previous_cancer == "no", "LFS-H",
                           ifelse(data_plot$cancer_status == "negative" & data_plot$previous_cancer == "yes", "LFS-PC",
                                  ifelse(data_plot$cancer_status == "positive", "LFS-AC", "Healthy")))
data_plot$diag <- factor(data_plot$diag, levels = c("HBC", "LFS"))
data_plot$status <- factor(data_plot$status, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))

### Calculate stats
data_stats_shan <- data_plot %>%
  group_by(status) %>%
  dplyr::summarise(mean=mean(shannon),
                   median=median(shannon),
                   sd=sd(shannon),
                   N=n())

stats <- unique(data_stats_shan$status)
x <- c()
for (stat in stats) {
  a <- t.test(data_plot$shannon[data_plot$status == "Healthy"], 
              data_plot$shannon[data_plot$status == stat])$p.value
  x <- c(x, a)
}
data_stats_shan$pvalue <- x
data_stats_shan$annot <- ifelse(data_stats_shan$pvalue < 0.05 & data_stats_shan$pvalue > 0.01, "*",
                                 ifelse(data_stats_shan$pvalue < 0.01 & data_stats_shan$pvalue > 0.001, "**",
                                        ifelse(data_stats_shan$pvalue < 0.001, "***", "")))

model <- aov(shannon ~ status, data = data_plot)
x <- summary(model)
y <- round(x[[1]][["Pr(>F)"]], 4)
anova_p_shannon <- y[[1]]

data_stats_gini <- data_plot %>%
  group_by(status) %>%
  dplyr::summarise(mean=mean(gini),
                   median=median(gini),
                   sd=sd(gini),
                   N=n())

stats <- unique(data_stats_gini$status)
x <- c()
for (stat in stats) {
  a <- t.test(data_plot$gini[data_plot$status == "Healthy"], 
              data_plot$gini[data_plot$status == stat])$p.value
  x <- c(x, a)
}
data_stats_gini$pvalue <- x
data_stats_gini$annot <- ifelse(data_stats_gini$pvalue < 0.05 & data_stats_gini$pvalue > 0.01, "*",
                                 ifelse(data_stats_gini$pvalue < 0.01 & data_stats_gini$pvalue > 0.001, "**",
                                        ifelse(data_stats_gini$pvalue < 0.001, "***", "")))

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               legend.position = "right",
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Plot comparisons
fig <- ggplot(data_plot) +
  geom_boxplot(aes(status, shannon, fill = status), alpha = 0.5, outlier.size = 1) +
  geom_text(data = data_stats_shan, aes(status, y = 0.9526, label = annot), size = 5) +
  geom_text(data = data_stats_shan, aes(status, y = 0.94, label = N), size = 4) +
  xlab("Status") + 
  ylab("Shannon") +
  labs(color = "", fill = "") +
  scale_fill_manual(name = " ", values =c("black", "#1F78B4", "#33A02C", "#E31A1C")) +
  ggtitle("Shannon Entropy") + 
  theme +
  scale_y_continuous(limits = c(0.939, 0.953))
fig

ggsave(file.path(outdir, paste0("fragment_end_contexts_shannon.pdf")), fig, width = 3, height = 4)

### Make Medians
lfsh <- data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "no"]
lfspc <- data_samples$sWGS[data_samples$cancer_status == "negative" & data_samples$previous_cancer == "yes"]
lfsac <- data_samples$sWGS[data_samples$cancer_status == "positive"]

hbc_median <- data.frame(motif = data_normal$motif,
                         median = 1,
                         type = "Healthy")
lfsh_median <- data.frame(motif = data$motif,
                          median = rowMedians(as.matrix(data[,colnames(data) %in% lfsh])/rowMedians(as.matrix(data_normal[,-1]))),
                          type = "LFS-H")
lfspc_median <- data.frame(motif = data$motif,
                           median = rowMedians(as.matrix(data[,colnames(data) %in% lfspc]))/rowMedians(as.matrix(data_normal[,-1])))),
                           type = "LFS-PC")
lfsac_median <- data.frame(motif = data$motif,
                           median = rowMedians(as.matrix(data[,colnames(data) %in% lfsac]/rowMedians(as.matrix(data_normal[,-1])))),
                           type = "LFS-AC")

data_plot_median <- bind_rows(hbc_median, lfsh_median, lfspc_median, lfsac_median)

my_chars <- strsplit(data_plot_median$motif, "")
char_frequency <- lapply(my_chars, function(x) {table(x)})
char_frequency <- do.call(bind_rows, char_frequency)
char_frequency <- as.data.frame(char_frequency)
char_frequency[is.na(char_frequency)] <- 0
char_frequency$AT <- char_frequency$A + char_frequency$T

data_plot_median$AT <- char_frequency$AT
data_plot_median$AT <- factor(data_plot_median$AT, levels = c(0,1,2,3,4))
order <- hbc_median$motif[order(hbc_median$median, decreasing = TRUE)]
data_plot_median$motif <- factor(data_plot_median$motif, levels = order)

### Plot medians
fig_median <- ggplot(data_plot_median) +
  geom_point(aes(motif, median, color = type), alpha = 0.5) +
  facet_grid(.~AT, scales = "free", space = "free") +
  xlab("Motif") + 
  ylab("Frequency") +
  labs(color = "", fill = "") +
  scale_color_manual(name = " ", values =c("black", "#1F78B4", "#33A02C", "#E31A1C")) +
  ggtitle("") + 
  theme +
  scale_y_log10()
fig_median

ggsave(file.path(outdir, paste0("fragment_end_contexts_shannon.pdf")), fig, width = 3, height = 4)
