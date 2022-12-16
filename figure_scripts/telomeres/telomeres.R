library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpubr)
library(lemon)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/telomeres"
source("//Users/derekwong/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Read in data
samples <- read.delim(file.path(path, "samples/sample_list.txt"))
data_contexts <- read.delim(list.files(path, "TelomereHunter_contexts_normalized.txt", recursive = TRUE, full.names = TRUE))
data_percents <- read.delim(list.files(path, "TelomereHunter_contexts_percent.txt", recursive = TRUE, full.names = TRUE))
data_summaries <- read.delim(list.files(path, "TelomereHunter_summaries.txt", recursive = TRUE, full.names = TRUE))
data_telseq <- read.delim(list.files(path, "TelSeq.txt", recursive = TRUE, full.names = TRUE))
data_prop <- read.delim(list.files(path, "proportions.txt", recursive = TRUE, full.names = TRUE))

healthy_samples <- read.delim(list.files(healthy_path, "sample_list.txt", recursive = TRUE, full.names = TRUE))
healthy_contexts <- read.delim(list.files(healthy_path, "TelomereHunter_contexts_normalized.txt", recursive = TRUE, full.names = TRUE))
healthy_percents <- read.delim(list.files(healthy_path, "TelomereHunter_contexts_percent.txt", recursive = TRUE, full.names = TRUE))
healthy_summaries <- read.delim(list.files(healthy_path, "TelomereHunter_summaries.txt", recursive = TRUE, full.names = TRUE))
healthy_telseq <- read.delim(list.files(healthy_path, "TelSeq.txt", recursive = TRUE, full.names = TRUE))
healthy_prop <- read.delim(list.files(healthy_path, "proportions.txt", recursive = TRUE, full.names = TRUE))

### Remove failed samples and format sample list
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]
samples <- samples[samples$sWGS %in% data_telseq$Library, ]

healthy_samples <- healthy_samples[healthy_samples$sWGS %in% healthy_telseq$Library, ]
samples <- bind_rows(samples, healthy_samples)
samples[is.na(samples)] <- "HBC"

### Format sequence contexts from TelomereHunter
data_contexts <- data_summaries[, colnames(data_summaries) %in% c("sample", "tel_content")]
data_contexts$diag <- "LFS"

healthy_contexts <- healthy_summaries[, colnames(healthy_summaries) %in% c("sample", "tel_content")]
healthy_contexts$diag <- "HBC"

contexts <- dplyr::bind_rows(data_contexts, healthy_contexts)
contexts <- merge(contexts, samples, by.x = "sample", by.y = "sWGS")

contexts$cancer_status <- factor(contexts$cancer_status, levels = c("HBC", "negative", "positive"),
                                 labels = c("Healthy", "LFS Negative", "LFS Positive"))

### Format summaries from TelomereHunter
data_summaries$diag <- "LFS"
healthy_summaries$diag <- "HBC"

summaries <- dplyr::bind_rows(data_summaries, healthy_summaries)
summaries <- merge(summaries, samples, by.x = "sample", by.y = "sWGS")

### Format TelSeq data
data_telseq$diag <- "LFS"
healthy_telseq$diag <- "HBC"

telseq <- dplyr::bind_rows(data_telseq, healthy_telseq)
telseq <- merge(telseq, samples, by.x = "Library", by.y = "sWGS")

telseq$cancer_status <- factor(telseq$cancer_status, levels = c("HBC", "negative", "positive"),
                               labels = c("Healthy", "LFS Negative", "LFS Positive"))

### Set theme and colors
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13),
               axis.line = element_line(colour = "black"),
               axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.background = element_blank(),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))

col <- c("grey65", "#FB9A99", "#FB9A99", "#FB9A99", "#FB9A99", "#FB9A99", "#FB9A99")

### Create comparison matrix for TelomereHunter to TelSeq
comparison_telseq <- telseq[, colnames(telseq) %in% c("Library", "tel_norm", "diag")]
comparison_summaries <- summaries[, colnames(summaries) %in% c("sample", "tel_content")]
comparison <- merge(comparison_telseq, comparison_summaries, by.x = "Library", by.y = "sample")
comparison <- merge(comparison, samples, by.x = "Library", by.y = "sWGS", all = TRUE)
rm(comparison_telseq, comparison_summaries)

plot_comparison <- ggplot(comparison, aes(tel_content, tel_norm, color = diag)) +
  geom_point(pch = 16, alpha = 1) +
  geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = 0.5, size = 1) +
  stat_regline_equation(aes(label = ..rr.label..), label.y = c(100, 50, 25), label.x = 100, size = 5) +
  scale_color_manual(values = c("red", "black")) +
  xlab("TelomereHunter (TRPM)") + 
  ylab("TelSeq (TRPM)") +
  labs(color = "Patient Type") +
  ggtitle("") + 
  theme +
  theme(legend.position = "bottom")
plot_comparison

ggsave(file.path(outdir, "telomere_sup_1.pdf"), plot_comparison, device = "pdf", width = 5, height = 5, units = "in")

### Plot each analysis vs age
contexts$years <- as.numeric(contexts$years)
plot_telhun <- ggplot(contexts[!(contexts$diag == "HBC"), ], aes(years, tel_content, color = cancer_status)) +
  geom_point(pch = 16) +
  geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = 0.25, size = 1) +
  stat_regline_equation(label.y = c(200, 225), label.x = 10, aes(label = ..rr.label..)) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Age (Years)") + 
  ylab("Telomere Content (TRPM)") +
  ggtitle("Telomere Hunter") + 
  theme +
  scale_y_continuous(limits = c(0, 250))
plot_telhun

telseq$years <- as.numeric(telseq$years)
plot_telseq <- ggplot(telseq[!(telseq$diag == "HBC"), ], aes(years, tel_norm, color = cancer_status)) +
  geom_point(pch = 16) +
  geom_line(stat = "smooth", method = "lm", se = FALSE, alpha = 0.25, size = 1) +
  stat_regline_equation(label.y = c(450, 500), label.x = 40, aes(label = ..rr.label..)) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Age (Years)") + 
  ylab("Telomere Content (TRPM)") +
  ggtitle("TelSeq") + 
  theme +
  scale_y_continuous(limits = c(0, 550))
plot_telseq

Figure <- ggpubr::ggarrange(plot_telhun, plot_telseq, nrow = 1, align = "h", widths = c(1,1))
Figure
ggsave(file.path(outdir, "telomere_vs_age.pdf"), Figure, device = "pdf", width = 8, height = 4, units = "in")

### Calculate stats for comparisons
# TelomereHunter tel_content
stats_conts <- contexts %>%
  group_by(cancer_status) %>% 
  dplyr::summarise(Median=median(tel_content, na.rm = TRUE),
                   Mean=mean(tel_content, na.rm = TRUE),
                   SD=sd(tel_content, na.rm = TRUE),
                   N=n())
a <- t.test(contexts$tel_content[contexts$cancer_status == "Healthy"], contexts$tel_content[contexts$cancer_status == "LFS Negative"])$p.value
b <- t.test(contexts$tel_content[contexts$cancer_status == "Healthy"], contexts$tel_content[contexts$cancer_status == "LFS Positive"])$p.value
t_test <- c(1, a, b)
stats_conts$pvalue <- t_test
stats_conts$annot <- ifelse(stats_conts$pvalue < 0.05 & stats_conts$pvalue > 0.01, "*",
                           ifelse(stats_conts$pvalue < 0.01 & stats_conts$pvalue > 0.001, "**",
                                  ifelse(stats_conts$pvalue < 0.001, "***", "")))

# Telseq tel_norm
stats_tels <- telseq %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(tel_norm, na.rm = TRUE),
                   Mean=mean(tel_norm, na.rm = TRUE),
                   SD=sd(tel_norm, na.rm = TRUE),
                   N=n())
a <- t.test(telseq$tel_norm[telseq$cancer_status == "Healthy"], telseq$tel_norm[telseq$cancer_status == "LFS Negative"])$p.value
b <- t.test(telseq$tel_norm[telseq$cancer_status == "Healthy"], telseq$tel_norm[telseq$cancer_status == "LFS Positive"])$p.value
t_test <- c(1, a, b)
stats_tels$pvalue <- t_test
stats_tels$annot <- ifelse(stats_tels$pvalue < 0.05 & stats_tels$pvalue > 0.01, "*",
                            ifelse(stats_tels$pvalue < 0.01 & stats_tels$pvalue > 0.001, "**",
                                   ifelse(stats_tels$pvalue < 0.001, "***", "")))


### Plot different analyses
plot_cont <- ggplot(contexts, aes(cancer_status, tel_content, fill = cancer_status)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_conts, aes(x = cancer_status, y = -5, label = N), size = 4) +
  geom_text(data = stats_conts, aes(x = cancer_status, y = 225, label = annot), size = 5) +
  scale_fill_manual(values = c("grey", "#a6cee3", "#fb9a99")) +
  xlab("") + 
  ylab("Telomere Reads per Million") +
  ggtitle("Telomere Content") + 
  theme +
  scale_y_continuous(limits=c(-20, 200), expand = c(0,0))
plot_cont

plot_tel <- ggplot(telseq, aes(cancer_status, tel_norm, fill = cancer_status)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_tels, aes(x = cancer_status, y = -10, label = N), size = 4) +
  geom_text(data = stats_tels, aes(x = cancer_status, y = 525, label = annot), size = 5) +
  scale_fill_manual(values = c("grey", "#a6cee3", "#fb9a99")) +
  xlab("") + 
  ylab("Telomere Content (TRPM)") +
  ggtitle("TelSeq") + 
  theme +
  scale_y_continuous(limits=c(-50, 500), expand = c(0,0))
plot_tel

ggsave(file.path(outdir, "telomere_diagnosis_telhun.pdf"), plot_cont, device = "pdf", width = 2, height = 5, units = "in")
ggsave(file.path(outdir, "telomere_diagnosis_telseq.pdf"), plot_tel, device = "pdf", width = 2, height = 5, units = "in")

### Calculate stats for tel_content sequences from TelomereHunter vs clinical
# cancer status
stats_status <- summaries %>%
  group_by(cancer_status)%>% 
  dplyr::summarise(Median=median(tel_content, na.rm = TRUE),
                   Mean=mean(tel_content, na.rm = TRUE),
                   SD=sd(tel_content, na.rm = TRUE),
                   N=n())
my_comparisons_status <- list(c("HBC", "negative"))

# age
stats_age <- summaries %>%
  group_by(Age)%>% 
  dplyr::summarise(Median=median(tel_content, na.rm = TRUE),
                   Mean=mean(tel_content, na.rm = TRUE),
                   SD=sd(tel_content, na.rm = TRUE),
                   N=n())
my_comparisons_age <- list(c("HBC", "pediatric"),
                           c("pediatric", "adult"))

# previous cancer
summaries$previous_cancer <- factor(summaries$previous_cancer, levels = c("HBC", "no", "yes"))
stats_previous <- summaries[summaries$previous_cancer %in% c("HBC", "no", "yes"), ] %>%
  group_by(previous_cancer)%>% 
  dplyr::summarise(Median=median(tel_content, na.rm = TRUE),
                   Mean=mean(tel_content, na.rm = TRUE),
                   SD=sd(tel_content, na.rm = TRUE),
                   N=n())
my_comparisons_previous <- list(c("HBC", "no"),
                                c("HBC", "yes"))

### Plot comparisons
plot_status <- ggplot(summaries, aes(cancer_status, tel_content, fill = cancer_status)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_status, aes(x = cancer_status, y = -5, label = N), size = 4) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Telomere Content (TRPM)") +
  ggtitle("Cancer Status") + 
  theme +
  scale_y_continuous(limits=c(-15, 300), expand = c(0,0))
plot_status

summaries$Age <- factor(summaries$Age, levels = c("HBC", "pediatric", "adult"))
plot_age <- ggplot(summaries, aes(Age, tel_content, fill = Age)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_age, aes(x = Age, y = -5, label = N), size = 4) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Telomere Content (TRPM)") +
  ggtitle("Age") + 
  theme +
  scale_y_continuous(limits=c(-15, 300), expand = c(0,0))
plot_age

plot_previous <- ggplot(summaries[summaries$previous_cancer %in% c("HBC", "no", "yes"), ], aes(previous_cancer, tel_content, fill = previous_cancer)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = stats_previous, aes(x = previous_cancer, y = -5, label = N), size = 4) +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Telomere Content (TRPM)") +
  ggtitle("Cancer History") + 
  theme +
  scale_y_continuous(limits=c(-15, 300), expand = c(0,0))
plot_previous

Figure <- ggpubr::ggarrange(plot_age, plot_status, plot_previous, nrow = 1, align = "h")
Figure
ggsave(file.path(outdir, "telomere_scores.pdf"), Figure, device = "pdf", width = 6, height = 4, units = "in")

### Create comparison matrix to compare to fragment proportions
data_prop <- data_prop[data_prop$sample %in% samples$sWGS, ]
data_prop <- rbind(data_prop, healthy_prop)

comparison <- merge(comparison, data_prop[,c(1,4)], by.x = "Library", by.y = "sample")
comparison$cancer_status <- factor(comparison$cancer_status, levels = c("HBC", "negative", "positive"),
                                   labels = c("Healthy", "LFS Negative", "LFS Positive"))

### Plot comparisons of telomeres vs fragment proportion
LFS_median <- median(data_prop$P.20_150.)

plot_tel_content <- ggplot(comparison[comparison$diag == "LFS", ], aes(tel_content, `P.20_150.`, color = cancer_status)) +
  geom_point(stroke = 0, pch = 16, size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1, alpha = 0.5) +
  stat_regline_equation(label.y = c(0.05, 0.1), label.x = 100, aes(label = ..rr.label..), size = 5) +
  geom_hline(yintercept = LFS_median, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("black", "red", "red")) +
  xlab("Telomere Content (TRPM)") + 
  ylab("Proportion under 150bp") +
  labs(color = "Cancer Status") +
  ggtitle("Fragment Proportion vs Telomere Content") + 
  theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 0.45), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0, 200))
plot_tel_content

ggsave(file.path(outdir, "telomere_1.pdf"), plot_tel_content, device = "pdf", width = 4.5, height = 5, units = "in")

### Plot telomeres across time
comparison2 <- comparison[comparison$diag == "LFS", ]
comparison2 <- comparison2[comparison2$ext_ID %in% comparison2$ext_ID[duplicated(comparison2$ext_ID)], ]
comparison2 <- comparison2 %>%
  group_by(ext_ID) %>%
  arrange(timepoint)

plot_patients <- ggplot(comparison2, aes(ext_ID, tel_content, group = ext_ID)) +
  geom_line(position = position_dodge2(0.5)) +
  geom_point(aes(color = cancer_status), 
             alpha = 0.75, position = position_dodge2(0.5), size = 1.5, pch = 16) +
  scale_color_manual(labels = c("Negative", "Positive"), values = c("black", "red")) +
  xlab("Patient") + 
  ylab("Telomeres (TRPM)") +
  labs(color = "Cancer Status") +
  ggtitle("") + 
  theme +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1))
plot_patients

ggsave(file.path(outdir, "telomere_patient.pdf"), plot_patients, device = "pdf", width = 12, height = 2.5, units = "in")

### Create matrix to compare telomere contexts
data_percents$diag <- "LFS"
data_percents$other <- rowSums2(as.matrix(data_percents[, !(colnames(data_percents) %in% c("sample", "TTAGGG", "TCAGGG", "TGAGGG", "diag", "sums"))]))
data_percents <- data_percents[, colnames(data_percents) %in% c("sample", "TTAGGG", "TCAGGG", "TGAGGG", "other", "diag")]
data_percents <- merge(data_percents, samples, by.x = "sample", by.y = "sWGS")

healthy_percents$other <- rowSums2(as.matrix(healthy_percents[, !(colnames(healthy_percents) %in% c("sample", "TTAGGG", "TCAGGG", "TGAGGG", "diag", "sums"))]))
healthy_percents <- healthy_percents[, colnames(healthy_percents) %in% c("sample", "TTAGGG", "TCAGGG", "TGAGGG", "other")]

percents <- dplyr::bind_rows(data_percents, healthy_percents)
percents[is.na(percents)] <- "HBC"

### Format contexts for comparison with clinical and plot
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13),
               axis.line = element_line(colour = "black"),
               axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
               axis.text.x = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "right",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.background = element_blank(),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13))
col <- c("#D95F02", "#7570B3", "#E6AB02", "#66A61E")

data_percents <- data_percents[, 1:6]
healthy_percents$diag <- "HBC"
percents <- bind_rows(data_percents, healthy_percents)
percents_melt <- reshape2::melt(percents, id = c("sample", "diag"))
percents_melt <- merge(percents_melt, samples, by.x = "sample", by.y = "sWGS", all = TRUE)
percents_melt[is.na(percents_melt)] <- "HBC"
percents_melt$cancer_status <- factor(percents_melt$cancer_status, levels = c("HBC", "negative", "positive"))
percents_melt$Age <- factor(percents_melt$Age, levels = c("HBC", "pediatric", "adult"))
percents_melt$variable <- factor(percents_melt$variable, levels = c("other", "TCAGGG", "TGAGGG", "TTAGGG"))
percents_melt <- percents_melt[order(percents_melt$cancer_status,
                                     percents_melt$variable,
                                     percents_melt$value), ]

order <- percents_melt[percents_melt$variable == "TTAGGG", 1]
percents_melt$sample <- factor(percents_melt$sample, levels = order)

plot_status <- ggplot(percents_melt) +
  geom_bar(aes(sample, value, fill = variable), position = "stack", stat = "identity") +
  facet_grid(.~cancer_status, scales = "free", space = "free") +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Proportion (%)") +
  labs(fill = "Context") +
  ggtitle("Cancer Status") + 
  theme +
  scale_y_continuous(limits=c(0, 1.01), expand = c(0,0))
plot_status

percents_melt <- percents_melt[order(percents_melt$Age,
                                     percents_melt$variable,
                                     percents_melt$value), ]
order <- percents_melt[percents_melt$variable == "TTAGGG", 1]
percents_melt$sample <- factor(percents_melt$sample, levels = order)

plot_age <- ggplot(percents_melt) +
  geom_bar(aes(sample, value, fill = variable), position = "stack", stat = "identity") +
  facet_grid(.~Age, scales = "free", space = "free") +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Proportion (%)") +
  labs(fill = "Context") +
  ggtitle("Age") + 
  theme +
  scale_y_continuous(limits=c(0, 1.01), expand = c(0,0))
plot_age

percents_years <- percents_melt[percents_melt$diag == "LFS", ]
percents_years <- percents_years[order(percents_years$years,
                                       percents_years$variable,
                                       percents_years$value), ]
order <- percents_years[percents_years$variable == "TTAGGG", 1]
percents_years$sample <- factor(percents_years$sample, levels = order)

plot_years <- ggplot(percents_years) +
  geom_bar(aes(sample, value, fill = variable), position = "stack", stat = "identity") +
  scale_fill_manual(values = col) +
  xlab("") + 
  ylab("Proportion (%)") +
  labs(fill = "Context") +
  ggtitle("Age (Years)") + 
  theme +
  scale_y_continuous(limits=c(0, 1.01), expand = c(0,0))
plot_years

Figure <- ggpubr::ggarrange(plot_status, plot_age, plot_years, nrow = 3)
Figure
ggsave(file.path(outdir, "telomere_sup_2.pdf"), Figure, device = "pdf", width = 16, height = 8, units = "in")

