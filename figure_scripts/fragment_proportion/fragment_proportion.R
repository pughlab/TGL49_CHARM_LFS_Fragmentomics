library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpmisc)
library(ggpubr)
library(GGally)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/insert_size"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_proportion"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/insert_size"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
source("//Users/derekwong/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Find files
proportion <- read.delim(list.files(path, "proportion", full.names = TRUE))
normal_prop <- read.delim(list.files(healthy_path, "proportion", full.names = TRUE))
data_samples <- read.delim(samples)

### Remove failed and unknown samples and format 
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
proportion <- proportion[!(proportion$cancer == "unknown"), ]
proportion <- proportion[!(proportion$sample %in% exclude), ]
proportion <- proportion[proportion$sample %in% data_samples$sWGS, ]
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[data_samples$sWGS %in% proportion$sample, ]

## Adjust and set names to data
data_proportion <- rbind(proportion, normal_prop)
colnames(data_proportion) <- c("sample", "cancer", "type", "p_20_150")
data_proportion$cancer <- factor(data_proportion$cancer,
                      levels = c("healthy_CHARM", "negative", "positive"),
                      labels = c("Healthy", "Negative", "Positive"))
data_proportion$type <- factor(data_proportion$type, levels = c("healthy", "LFS"),
                               labels = c("Healthy", "LFS"))
data_all <- data_proportion[complete.cases(data_proportion), ]
data_prop <- data_proportion[complete.cases(data_proportion) & data_proportion$type %in% c("Healthy", "LFS"), ]
data_prop$cancer <- factor(data_prop$cancer, levels = c("Healthy", "Negative", "Positive"),
                           labels = c("Healthy", "LFS Negative", "LFS Positive"))

## Calculate statistics (All LFS)
data_stats_LFS <- data_prop %>%
  group_by(type) %>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
a <- t.test(data_prop$p_20_150[data_prop$type=="Healthy"], data_prop$p_20_150[data_prop$type=="LFS"])$p.value
t_test <- c(1, a)
data_stats_LFS$pvalue <- t_test
data_stats_LFS$annot <- ifelse(data_stats_LFS$pvalue < 0.05 & data_stats_LFS$pvalue > 0.01, "*",
                                  ifelse(data_stats_LFS$pvalue < 0.01 & data_stats_LFS$pvalue > 0.001, "**",
                                         ifelse(data_stats_LFS$pvalue < 0.001, "***", "")))

### Calculate statistics (cancer status)
data_stats_cancer <- data_prop %>%
  group_by(cancer)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
a <- t.test(data_prop$p_20_150[data_prop$cancer=="Healthy"], data_prop$p_20_150[data_prop$cancer=="LFS Negative"])$p.value
b <- t.test(data_prop$p_20_150[data_prop$cancer=="Healthy"], data_prop$p_20_150[data_prop$cancer=="LFS Positive"])$p.value
c <- t.test(data_prop$p_20_150[data_prop$cancer=="LFS Negative"], data_prop$p_20_150[data_prop$cancer=="LFS Positive"])$p.value
t_test <- c(1, a, b)
t_test2 <- c(1, 1, c)
data_stats_cancer$pvalue <- t_test
data_stats_cancer$pvalue2 <- t_test2
data_stats_cancer$annot <- ifelse(data_stats_cancer$pvalue < 0.05 & data_stats_cancer$pvalue > 0.01, "*",
                           ifelse(data_stats_cancer$pvalue < 0.01 & data_stats_cancer$pvalue > 0.001, "**",
                                  ifelse(data_stats_cancer$pvalue < 0.001, "***", "")))
data_stats_cancer$annot2 <- ifelse(data_stats_cancer$pvalue2 < 0.05 & data_stats_cancer$pvalue2 > 0.01, "*",
                                  ifelse(data_stats_cancer$pvalue2 < 0.01 & data_stats_cancer$pvalue2 > 0.001, "**",
                                         ifelse(data_stats_cancer$pvalue2 < 0.001, "***", "")))
rm(a,b,c,t_test,t_test2)

my_comparisons <- list( c("Healthy", "LFS Negative"), 
                        c("Healthy", "LFS Positive"),
                        c("LFS Negative", "LFS Positive"))

### Calculate statistics (previvors)
samples_neg <- data_samples[data_samples$cancer_status == "negative" & 
                              data_samples$previous_cancer == "no", ]
data_prop_hbc <- data_prop[data_prop$type == "Healthy", ]
data_prop_neg <- data_prop[data_prop$sample %in% samples_neg$sWGS, ]
data_prop_neg <- bind_rows(data_prop_neg, data_prop_hbc)
data_prop_neg$type <- factor(data_prop_neg$type, levels = c("Healthy", "LFS"),
                             labels = c("Healthy", "Previvor"))

data_stats_neg <- data_prop_neg %>%
  group_by(type)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
a <- t.test(data_prop_neg$p_20_150[data_prop_neg$type =="Healthy"], data_prop_neg$p_20_150[data_prop_neg$type =="Previvor"])$p.value
t_test <- c(1, a)
data_stats_neg$pvalue <- t_test
data_stats_neg$annot <- ifelse(data_stats_neg$pvalue < 0.05 & data_stats_neg$pvalue > 0.01, "*",
                               ifelse(data_stats_neg$pvalue < 0.01 & data_stats_neg$pvalue > 0.001, "**",
                                      ifelse(data_stats_neg$pvalue < 0.001, "***", "")))

### Calculate statistics (cancer history)
data_prop_LFS <- data_prop[data_prop$type == "LFS", ]
data_prop_LFS <- merge(data_prop_LFS, data_samples, by.x = "sample", by.y = "sWGS")
data_stats_previous <- data_prop_LFS[!(data_prop_LFS$previous_cancer == ""), ] %>%
  group_by(previous_cancer)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
a <- t.test(data_prop_LFS$p_20_150[data_prop_LFS$previous_cancer=="no"], data_prop_LFS$p_20_150[data_prop_LFS$previous_cancer=="yes"])$p.value
t_test <- c(1, a)
data_stats_previous$pvalue <- t_test

### Calculate statistics (mutation type)
data_prop_type <- data_prop_LFS[data_prop_LFS$mutation_type %in% c("LOF", "1", "2", "3", "4", "5", "Splice"), ]
data_prop_type$mutation_type <- factor(data_prop_type$mutation_type, levels = c("1", "2", "3", "5", "LOF", "Splice"),
                                       labels = c("Missense 1", "Missense 2", "Missense 3", "Missense 5", "LOF", "Splice"))
data_stats_type <- data_prop_type %>%
  group_by(mutation_type)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n(),
                   patients=length(unique(sample_parent)))
t_test <- c()
vars <- data_stats_type$mutation_type
for (var in vars) {
  a <- t.test(data_prop$p_20_150[data_prop$type == "LFS"], data_prop_type$p_20_150[data_prop_type$mutation_type == var])$p.value
  t_test <- c(t_test, a)
}

data_stats_type$pvalue <- t_test
data_stats_type$annot <- ifelse(data_stats_type$pvalue < 0.05 & data_stats_type$pvalue > 0.01, "*",
                               ifelse(data_stats_type$pvalue < 0.01 & data_stats_type$pvalue > 0.001, "**",
                                      ifelse(data_stats_type$pvalue < 0.001, "***", "")))

### Calculate statistics (Germline Mutations)
data_prop_germline <- data_prop_LFS
data_prop_germline <- data_prop_germline[!(data_prop_germline$germline_mutation == "unknown"), ]
data_prop_germline$germline_mutation<- ifelse(data_prop_germline$germline_mutation %in% c("Exon 1 del", "intron8_deletion"), "Deletion", 
                                              data_prop_germline$germline_mutation )
data_prop_germline$germline_mutation<- ifelse(data_prop_germline$germline_mutation == "C275* A276_C277delinsG", "C275*", 
                                              data_prop_germline$germline_mutation )
counts <- as.data.frame(table(data_prop_germline$germline_mutation))
singles <- counts[counts$Freq < 2, ]
singles <- singles$Var1
data_prop_germline <- data_prop_germline[!(data_prop_germline$germline_mutation %in% singles), ]
data_prop_germline$mutation_type <- factor(data_prop_germline$mutation_type, levels = c("1", "2", "3", "5", "LOF", "Splice"),
                                           labels = c("Missense 1", "Missense 2", "Missense 3", "Missense 5", "LOF", "Splice"))

data_stats_germ <- data_prop_germline %>%
  group_by(germline_mutation)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n(),
                   patients=length(unique(sample_parent)),
                   type=unique(mutation_type))
t_test <- c()
vars <- data_stats_germ$germline_mutation
for (var in vars) {
  a <- t.test(data_prop$p_20_150[data_prop$type == "LFS"], data_prop_germline$p_20_150[data_prop_germline$germline_mutation == var])$p.value
  t_test <- c(t_test, a)
}

data_stats_germ$pvalue <- t_test
data_stats_germ$annot <- ifelse(data_stats_germ$pvalue < 0.05 & data_stats_germ$pvalue > 0.01, "*",
                                ifelse(data_stats_germ$pvalue < 0.01 & data_stats_germ$pvalue > 0.001, "**",
                                       ifelse(data_stats_germ$pvalue < 0.001, "***", "")))

data_stats_germ <- data_stats_germ[order(data_stats_germ$type,
                                         data_stats_germ$Median), ]
order <- data_stats_germ$germline_mutation
data_prop_germline$germline_mutation <- factor(data_prop_germline$germline_mutation, levels = order)

### Calculate statistics (families)
data_prop_fam <- data_prop_LFS
data_prop_fam$family <- gsub("LIB-04-", "", data_prop_fam$family)
data_prop_fam$family <- ifelse(data_prop_fam$family %in% c("5946", "5952", ""), "Remainder", data_prop_fam$family)
data_prop_fam$family <- ifelse(!(data_prop_fam$family == "Remainder"), 
                               paste0(data_prop_fam$family, " (", data_prop_fam$germline_mutation, ")"),
                               data_prop_fam$family)
data_prop_fam$germline_mutation <- ifelse(data_prop_fam$family == "Remainder", "Mixed", data_prop_fam$germline_mutation)
data_prop_fam$mutation_type <- ifelse(data_prop_fam$germline_mutation == "Mixed", "Mixed", data_prop_fam$mutation_type)
data_prop_fam$mutation_type <- factor(data_prop_fam$mutation_type, levels = c("1", "2", "3", "5", "LOF", "Splice", "Mixed"),
                                      labels = c("Missense 1", "Missense 2", "Missense 3", "Missense 5", "LOF", "Splice", "Mixed"))
data_stats_fam <- data_prop_fam %>%
  group_by(family)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n(),
                   patients=length(unique(sample_parent)),
                   mutation=unique(germline_mutation),
                   type=unique(mutation_type))

t_test <- c()
vars <- data_stats_fam$family
for (var in vars) {
  a <- t.test(data_prop$p_20_150[data_prop$type == "LFS"], data_prop_fam$p_20_150[data_prop_fam$family == var])$p.value
  t_test <- c(t_test, a)
}

data_stats_fam$pvalue <- t_test
data_stats_fam$annot <- ifelse(data_stats_fam$pvalue < 0.05 & data_stats_fam$pvalue > 0.01, "*",
                               ifelse(data_stats_fam$pvalue < 0.01 & data_stats_fam$pvalue > 0.001, "**",
                                      ifelse(data_stats_fam$pvalue < 0.001, "***", "")))

duplicated <- unique(data_stats_fam$mutation[duplicated(data_stats_fam$mutation)])
data_stats_fam$text <- ifelse(data_stats_fam$mutation %in% duplicated, "bold", "plain")

data_stats_fam <- data_stats_fam[order(data_stats_fam$type,
                                       data_stats_fam$Median), ]
order <- data_stats_fam$family
data_prop_fam$family <- factor(data_prop_fam$family, levels = order)

### Make healthy medians
healthy_median <- as.double(data_stats_cancer[data_stats_cancer$cancer == "Healthy", colnames(data_stats_cancer) == "Median"])
previvor_median <- as.double(data_stats_neg[data_stats_neg$type == "Previvor", colnames(data_stats_neg) == "Median"])
LFS_median <- as.double(data_stats_LFS[data_stats_LFS$type == "LFS", colnames(data_stats_neg) == "Median"])

## Plot Figures
### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 12),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Healthy vs LFS
fig <- ggplot(data_prop, aes(type, p_20_150) ) +
  geom_boxplot(aes(fill = type), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = type), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_hline(yintercept = healthy_median, linetype = "dashed", color = "red") +
  geom_text(data = data_stats_LFS, aes(x = type, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_LFS, aes(x = type, y = 0.5, label = annot), size = 5.5) +
  scale_fill_manual(values = c("grey", "#fb9a99")) +
  labs(fill = "") +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("") + 
  theme +
  scale_y_continuous(limits=c(0, 0.55), expand = c(0,0))
fig  

### Previvors
fig_neg <- ggplot(data_prop_neg, aes(type, p_20_150) ) +
  geom_boxplot(aes(fill = type), outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(aes(fill = type), color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_hline(yintercept = healthy_median, linetype = "dashed", color = "red") +
  geom_text(data = data_stats_neg, aes(x = type, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_neg, aes(x = type, y = 0.5, label = annot), size = 5.5) +
  scale_fill_manual(values = c("grey", "#fb9a99")) +
  labs(fill = "") +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("") + 
  theme +
  scale_y_continuous(limits=c(0, 0.55), expand = c(0,0))
fig_neg  

### Within LFS
fig_lfs <- ggplot(data_prop, aes(cancer, p_20_150, fill = cancer)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_hline(yintercept = healthy_median, linetype = "dashed", color = "red") +
  geom_text(data = data_stats_cancer, aes(x = cancer, y = 0.03, label = N), size = 4) +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("") + 
  scale_fill_manual(labels = c("Healthy", "Negative", "Positive"), 
                    values = c("grey", "#a6cee3", "#fb9a99")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(0, 0.55), expand = c(0,0)) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     size = 5,
                     tip.length = 0,
                     step.increase = 0.11,
                     position = 0.2)
fig_lfs

### Previous Cancer
fig_previous <- ggplot(data_prop_LFS[!(data_prop_LFS$previous_cancer == ""), ], aes(previous_cancer, p_20_150, fill = previous_cancer)) +
  geom_boxplot(outlier.size = 0.5, 
               position = position_nudge(x = 0.11), width = 0.4) +
  geom_flat_violin(color = NA, position = position_nudge(x = -0.11), width = 0.5) +
  geom_text(data = data_stats_previous, aes(x = previous_cancer, y = 0.03, label = N), size = 4) +
  geom_hline(yintercept = healthy_median, linetype = "dashed", color = "red") +
  xlab("Previous\nCancer") + 
  ylab("Proportion under 150bp") +
  ggtitle("") + 
  scale_fill_manual(values = c("grey", "#fb9a99")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(0, 0.55), expand = c(0,0))
fig_previous

### Age
data_prop_LFS$years <- as.numeric(data_prop_LFS$years)

fig_age <- ggplot(data_prop_LFS, aes(years, p_20_150, color = cancer_status)) + 
  geom_point(stroke = 0, pch = 16, size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1, alpha = 0.5) +
  stat_regline_equation(label.y = c(0.05, 0.1), label.x = 10, aes(label = ..rr.label..), size = 5) +
  geom_hline(yintercept = LFS_median, linetype = "dashed", color = "red") +
  ggtitle("Fragment vs Age") +
  xlab("Age (years)") +
  ylab("Proportion under 150bp") +
  labs(fill = "Cancer Status") +
  theme(plot.title = element_text(hjust = 0.5, size = 13), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.75, 0.2),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(limits = c(0, 0.4), expand = c(0,0)) + 
  scale_x_continuous(limits = c(0, 70)) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE))
fig_age

### Mutation Type
data_prop_type <- data_prop_type %>%
  group_by(mutation_type) %>%
  arrange(p_20_150)
data_stats_type$Xn <- c(1:nrow(data_stats_type))
fig_type <- ggplot(data_prop_type) +
  stat_summary(aes(mutation_type, p_20_150), 
               fun = median, geom = "crossbar", size = 0.2, width = 0.5, group = 1, show.legend = FALSE, color = "red") +
  geom_rect(data = data_stats_type, aes(xmin = Xn - 0.5, xmax = Xn + 0.5, ymin = -Inf, ymax = Inf, fill = mutation_type), alpha = 0.5, stat = "identity") +
  geom_point(aes(mutation_type, p_20_150, group = mutation_type, color = cancer), alpha = 0.75,
             position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_type, aes(x = mutation_type, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_type, aes(x = mutation_type, y = 0.06, label = paste0("(", patients, ")")), size = 4) +
  geom_text(data = data_stats_type, aes(x = mutation_type, y = 0.41, label = annot), size = 5) +
  geom_hline(yintercept = healthy_median, linetype = "dashed", color = "red", size = 0.25) +
  ggtitle("Germline Mutation Type") + 
  xlab("") +
  ylab("Proportion under 150bp") +
  scale_fill_manual(values = c("#FB9A99", "#A6CEE3", "#FDBF6F","#CAB2D6", "grey65", "#B2DF8A")) +
  scale_color_manual(values = c("black", "red")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(0, 0.45), expand = c(0,0))
fig_type

### Germline mutation
data_prop_germline <- data_prop_germline %>%
  group_by(germline_mutation) %>%
  arrange(p_20_150)
data_stats_germ$Xn <- c(1:nrow(data_stats_germ))
fig_germ <- ggplot(data_prop_germline) +
  stat_summary(aes(germline_mutation, p_20_150), 
               fun = median, geom = "crossbar", size = 0.2, width = 0.5, group = 1, show.legend = FALSE, color = "red") +
  geom_rect(data = data_stats_germ, aes(xmin = Xn - 0.5, xmax = Xn + 0.5, ymin = -Inf, ymax = Inf, fill = type), alpha = 0.5, stat = "identity") +
  geom_point(aes(germline_mutation, p_20_150, group = germline_mutation, color = cancer), alpha = 0.75,
             position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_germ, aes(x = germline_mutation, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_germ, aes(x = germline_mutation, y = 0.06, label = paste0("(", patients, ")")), size = 4) +
  geom_text(data = data_stats_germ, aes(x = germline_mutation, y = 0.41, label = annot), size = 5) +
  geom_hline(yintercept = healthy_median, linetype = "dashed", color = "red", size = 0.25) +
  ggtitle("Germline Mutation") + 
  xlab("") +
  ylab("Proportion under 150bp") +
  scale_fill_manual(values = c("#FB9A99", "#A6CEE3", "#FDBF6F", "#CAB2D6", "grey65", "#B2DF8A")) +
  scale_color_manual(values = c("black", "red")) +
  labs(fill = "") +
  theme +
  scale_y_continuous(limits=c(0, 0.45), expand = c(0,0))
fig_germ

### Family
data_prop_fam <- data_prop_fam %>%
  group_by(family) %>%
  arrange(p_20_150)
data_stats_fam$Xn <- c(1:nrow(data_stats_fam))
fig_fam <- ggplot(data_prop_fam) +
  stat_summary(aes(family, p_20_150), 
               fun = median, geom = "crossbar", size = 0.2, width = 0.5, group = 1, show.legend = FALSE, color = "red") +
  geom_rect(data = data_stats_fam, aes(xmin = Xn - 0.5, xmax = Xn + 0.5, ymin = -Inf, ymax = Inf, fill = type), alpha = 0.5, stat = "identity") +
  geom_point(aes(family, p_20_150, group = family, color = cancer), alpha = 0.75,
             position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_fam, aes(x = family, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_fam, aes(x = family, y = 0.06, label = paste0("(", patients, ")")), size = 4) +
  geom_text(data = data_stats_fam, aes(x = family, y = 0.41, label = annot), size = 5) +
  geom_hline(yintercept = healthy_median, linetype = "dashed", color = "red", size = 0.25) +
  ggtitle("LFS Family") + 
  xlab("") +
  ylab("Proportion under 150bp") +
  scale_fill_manual(breaks = c("Missense 1", "Missense 2", "Missense 3", "Missense 5", "LOF", "Splice", "Mixed"),
                    values = c("#FB9A99", "#A6CEE3", "#FDBF6F","#CAB2D6", "grey65", "#B2DF8A", "white")) +
  scale_color_manual(labels = c("Negative", "Positive"),
                     values = c("black", "red")) +
  labs(fill = "Mutation Type", color = "Cancer Status") +
  theme +
  theme(legend.position = "right",
        axis.text.x = element_text(face = data_stats_fam$text)) +
  scale_y_continuous(limits=c(0, 0.45), expand = c(0,0))
fig_fam

ggsave(file.path(outdir, "fragment_proportions_2.pdf"), fig_age, device = "pdf", width = 4.5, height = 5, units = "in")

ggsave(file.path(outdir, "fragment_proportions_3.pdf"), fig_type, device = "pdf", width = 4.5, height = 5, units = "in")

Figure <- ggarrange(fig_germ, fig_fam, align = "h", widths = c(1, 1), nrow = 1)
Figure
ggsave(file.path(outdir, "fragment_proportions_4.pdf"), Figure, device = "pdf", width = 16, height = 5, units = "in")

Figure_sup <- ggarrange(fig, fig_previous, align = "h", widths = c(1, 1, 3), nrow = 1)
Figure_sup
ggsave(file.path(outdir, "fragment_proportions_sup.pdf"), Figure_sup, device = "pdf", width = 4.5, height = 5, units = "in")


