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

### Find paths
data <- list.files(path, "dinucleotide.txt", full.names = TRUE)
data <- data[grepl("short", data)]
data_normal <- list.files(healthy_path, "dinucleotide.txt", full.names = TRUE)
data_normal <- data_normal[grepl("short", data_normal)]

### Import data 
data <- read.delim(data)
data_normal <- read.delim(data_normal)
data_samples <- read.delim("sample_list.txt")

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
#abnormal <- c("TGL49_0019_Cf_U_PE_326_WG", "TGL49_0022_Cf_U_PE_316_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% c(exclude)), ]
data_samples <- data_samples[data_samples$sWGS %in% colnames(data), ]
data <- data[, colnames(data) %in% c("context", data_samples$sWGS)]

### Sum frequencies across A/T and G/C combinations
data$context <- ifelse(data$context %in% c("AA", "AT", "TA", "TT"), "at",
                       ifelse(data$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
data <- data[complete.cases(data), ]
data_contexts <- aggregate(.~context, data, sum)

data_normal$context <- ifelse(data_normal$context %in% c("AA", "AT", "TA", "TT"), "at",
                              ifelse(data_normal$context %in% c("CC", "CG", "GC", "GG"), "gc", NA))
data_normal <- data_normal[complete.cases(data_normal), ]
data_normal_contexts <- aggregate(.~context, data_normal, sum)

### Make plotting table
data_melt <- reshape2::melt(data_contexts, id = "context")
data_melt <- merge(data_melt, data_samples, by.x = "variable", by.y = "sWGS")
data_melt$diag <- ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "no", "LFS-H",
                         ifelse(data_melt$cancer_status == "negative" & data_melt$previous_cancer == "yes", "LFS-PC", "LFS-AC"))

normal_melt <- reshape2::melt(data_normal_contexts, id = "context")

data <- bind_rows(data_melt, normal_melt)
data[is.na(data)] <- "Healthy"

### Format plotting table
data$diag <- factor(data$diag, levels = c("Healthy", "LFS-H", "LFS-PC", "LFS-AC"))
data$context <- factor(data$context, levels = c("at", "gc"),
                       labels = c("AA/AT/TA/TT", "CC/CG/GC/GG"))

### Calculate statistics
data_stats <- data %>%
  group_by(context, diag) %>%
  dplyr::summarise(mean=mean(value),
                   sd=sd(value),
                   N=n())
t_test <- c()
contexts <- unique(data_stats$context)
diags <- unique(data_stats$diag)
for (context in contexts) {
  for (diag in diags) {
    a <- t.test(data$value[data$context == context & data$diag == "Healthy"], data$value[data$context == context & data$diag == diag])$p.value
    t_test <- c(t_test, a)
  }
}
data_stats$pvalue <- t_test
data_stats$annot <- ifelse(data_stats$pvalue < 0.05 & data_stats$pvalue > 0.01, "*",
                           ifelse(data_stats$pvalue < 0.01 & data_stats$pvalue > 0.001, "**",
                                  ifelse(data_stats$pvalue < 0.001, "***", "")))

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill = NA),
               panel.background = element_blank(),
               legend.position = c(0.2, 0.35),
               legend.key = element_rect(fill = "white"),
               legend.title = element_text(size = 12),
               legend.text = element_text(size = 12),
               strip.background = element_blank(),
               strip.text = element_text(size = 13),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Plot curves
fig <- ggplot(data) +
  geom_boxplot(aes(diag, value, fill = diag), alpha = 0.5, outlier.size = 0.5) +
  geom_text(data = data_stats, aes(diag, 0.35, label = annot), size = 5) +
  geom_text(data = data_stats, aes(diag, 0.13, label = N), size = 4) +
  xlab("") + 
  ylab("Frequency") +
  labs(fill = "") +
  scale_fill_manual(values = c("black", "#1F78B4", "#33A02C", "#E31A1C")) +
  ggtitle("Fragments <150bp") + 
  facet_grid(.~context, scales = "free_y") +
  theme
fig

ggsave(file.path(outdir, "dinucleotide_contexts_short.pdf"), fig, width = 4, height = 5)

