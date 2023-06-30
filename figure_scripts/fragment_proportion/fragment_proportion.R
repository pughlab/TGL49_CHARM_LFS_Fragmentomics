library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpmisc)
library(ggpubr)
library(GGally)
library(patchwork)
library(lubridate)
library(data.table)

### Set paths
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/insert_size"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_proportion"
healthy_path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/insert_size"
samples <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
source("//Users/derekwong/My Drive/Post-Doc/CHARM/LFS/Figures/TP53_griffin/geom_flat_violin.R")

### Find files
proportion <- read.delim(list.files(path, "proportion", full.names = TRUE))
normal_prop <- read.delim(list.files(healthy_path, "proportion", full.names = TRUE))
groups <- read.delim(list.files("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS/read_group", ".txt", full.names = TRUE))
groups_n <- read.delim(list.files("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/read_group", ".txt", full.names = TRUE))
clusters <- read.delim("/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_proportion/cluster_density.txt")
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

model <- aov(p_20_150 ~ cancer, data = data_prop)
anova_result <- summary(model)
anova_p_status <- round(anova_result[[1]][["Pr(>F)"]], 4)
anova_p_status <- anova_p_status[[1]]

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

model <- aov(p_20_150 ~ mutation_type, data = data_prop_LFS)
anova_result <- summary(model)
anova_p_type <- round(anova_result[[1]][["Pr(>F)"]], 4)
anova_p_type <- anova_p_type[[1]]

model <- aov(p_20_150 ~ mutation_type, data = data_prop_LFS[data_prop_LFS$mutation_type %in% c("2", "3", "5", "Splice", "LOF"), ])
anova_result <- summary(model)
anova_p_type_minus <- round(anova_result[[1]][["Pr(>F)"]], 4)
anova_p_type_minus <- anova_p_type_minus[[1]]

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
  a <- wilcox.test(data_prop$p_20_150[data_prop$type == "LFS"], data_prop_germline$p_20_150[data_prop_germline$germline_mutation == var])$p.value
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

anova <- c()
kruskal <- c()
types <- c("2", "3", "Splice", "LOF")
for (type in types) {
  model <- aov(p_20_150 ~ germline_mutation, data = data_prop_LFS[data_prop_LFS$mutation_type == type, ])
  x <- summary(model)
  y <- round(x[[1]][["Pr(>F)"]], 4)
  y <- y[[1]]
  anova <- c(anova, y)
  
  model <- kruskal.test(p_20_150 ~ germline_mutation, data = data_prop_LFS[data_prop_LFS$mutation_type == type, ])
  a <- round(model$p.value, 4)
  kruskal <- c(kruskal, a)
}
anova_p_germ <- data.frame(germline_mutation = types,
                           pvalue = anova,
                           kruskal = kruskal)
model <- aov(p_20_150 ~ germline_mutation, data = data_prop_LFS)
x <- summary(model)
y <- round(x[[1]][["Pr(>F)"]], 4)
y <- y[[1]]

model <- kruskal.test(p_20_150 ~ germline_mutation, data = data_prop_LFS)
a <- round(model$p.value, 4)
anova_p_germ <- rbind(anova_p_germ, c("All", y, a))

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
  a <- wilcox.test(data_prop$p_20_150[data_prop$type == "LFS"], data_prop_fam$p_20_150[data_prop_fam$family == var])$p.value
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

model <- kruskal.test(p_20_150 ~ family, data = data_prop_LFS)
anova_p_fam <- round(model$p.value, 4)

data_R213Q <- data_prop_LFS[data_prop_LFS$germline_mutation == "R213Q", ]
model <- kruskal.test(p_20_150 ~ family, data = data_R213Q)
anova_p_R213Q <- round(model$p.value, 4)

data_T125 <- data_prop_LFS[data_prop_LFS$germline_mutation == "T125=_splice", ]
model <- kruskal.test(p_20_150 ~ family, data = data_T125)
anova_p_T125 <- round(model$p.value, 4)

### Calculate stats based on institution
data_prop_LFS$institution <- ifelse(data_prop_LFS$LIB_ID %like% "LFS", "SK", "PM")

data_stats_inst <- data_prop_LFS %>%
  group_by(institution)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
data_stats_inst$pvalue = c(1, t.test(data_prop_LFS$p_20_150[data_prop_LFS$institution == "SK"], data_prop_LFS$p_20_150[data_prop_LFS$institution == "PM"])$p.value)

### Make healthy medians
healthy_median <- as.double(data_stats_cancer[data_stats_cancer$cancer == "Healthy", colnames(data_stats_cancer) == "Median"])
previvor_median <- as.double(data_stats_neg[data_stats_neg$type == "Previvor", colnames(data_stats_neg) == "Median"])
LFS_median <- as.double(data_stats_LFS[data_stats_LFS$type == "LFS", colnames(data_stats_neg) == "Median"])
cohort_median <- median(data_prop$p_20_150)

### Format cluster information
clusters$run.ID <- paste0("PU:", clusters$run.ID)
clusters <- clusters[clusters$run.ID %in% c(groups$V1, groups_n$V1), ]

### Attach read groups
groups <- groups[groups$sample %in% data_prop$sample, ]
groups <- bind_rows(groups, groups_n)
data_reads <- merge(data_prop, groups, by.x = "sample", by.y = "sample", all = TRUE)
read_groups <- unique(data_reads$V1)
data_reads$V1 <- factor(data_reads$V1, levels = read_groups,
                        labels = paste0("run", c(1:length(read_groups))))

clusters$run.ID <- factor(clusters$run.ID, levels = read_groups,
                          labels = paste0("run", c(1:length(read_groups))))

### Order by cluster density
cluster_order <- clusters$run.ID[order(clusters$X, decreasing = TRUE)]
clusters$run.ID <- factor(clusters$run.ID, levels = cluster_order)

data_reads$V1 <- factor(data_reads$V1, levels = cluster_order)

### Calculate sequencing stats
data_stats_group <- data_reads %>%
  group_by(V1)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
                   Mean=mean(p_20_150, na.rm = TRUE),
                   SD=sd(p_20_150, na.rm = TRUE),
                   N=n())
data_stats_group <- data_stats_group[data_stats_group$N > 2, ]
data_reads <- data_reads[data_reads$V1 %in% data_stats_group$V1, ]
clusters <- clusters[clusters$run.ID %in% data_stats_group$V1, ]

t_test <- c()
vars <- data_stats_group$V1
for (var in vars) {
  a <- t.test(data_prop$p_20_150, data_reads$p_20_150[data_reads$V1 == var])$p.value
  t_test <- c(t_test, a)
}

data_stats_group$pvalue <- t_test
data_stats_group$pvalue <- p.adjust(data_stats_group$pvalue)
data_stats_group$annot <- ifelse(data_stats_group$pvalue < 0.05 & data_stats_group$pvalue > 0.01, "*",
                               ifelse(data_stats_group$pvalue < 0.01 & data_stats_group$pvalue > 0.001, "**",
                                      ifelse(data_stats_group$pvalue < 0.001, "***", "")))

model <- aov(p_20_150 ~ V1, data = data_reads)
anova_result <- summary(model)
anova_p <- round(anova_result[[1]][["Pr(>F)"]], 4)
anova_p <- anova_p[[1]]

data_reads <- merge(data_reads, clusters, by.x = "V1", by.y = "run.ID", all = TRUE)

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
  ggtitle("All LFS") + 
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
  ggtitle("LFS-H") + 
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
  ggtitle("LFS Cancer Status") + 
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
  ggtitle("Cancer History") + 
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
  #geom_rect(data = data_stats_type, aes(xmin = Xn - 0.5, xmax = Xn + 0.5, ymin = -Inf, ymax = Inf, fill = mutation_type), alpha = 0.5, stat = "identity") +
  geom_boxplot(aes(mutation_type, p_20_150, fill = mutation_type), alpha = 0.5, outlier.shape = NA) +
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
  #geom_rect(data = data_stats_germ, aes(xmin = Xn - 0.5, xmax = Xn + 0.5, ymin = -Inf, ymax = Inf, fill = type), alpha = 0.5, stat = "identity") +
  geom_boxplot(aes(germline_mutation, p_20_150, fill = mutation_type), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(germline_mutation, p_20_150, group = germline_mutation, color = cancer), alpha = 0.75,
             position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_germ, aes(x = germline_mutation, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_germ, aes(x = germline_mutation, y = 0.06, label = paste0("(", patients, ")")), size = 4) +
  geom_text(data = data_stats_germ, aes(x = germline_mutation, y = 0.41, label = annot), size = 5) +
  geom_hline(yintercept = LFS_median, linetype = "dashed", color = "red", size = 0.25) +
  ggtitle("Germline Mutation") + 
  xlab("") +
  ylab("Proportion under 150bp") +
  scale_fill_manual(values = c("#FB9A99", "#A6CEE3", "#FDBF6F", "#CAB2D6", "grey65", "#B2DF8A")) +
  scale_color_manual(values = c("black", "red")) +
  labs(fill = "Mutation Type", color = "Cancer Status") +
  theme +
  theme(legend.position = "right") +
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
  #geom_rect(data = data_stats_fam, aes(xmin = Xn - 0.5, xmax = Xn + 0.5, ymin = -Inf, ymax = Inf, fill = type), alpha = 0.5, stat = "identity") +
  geom_boxplot(aes(family, p_20_150, fill = mutation_type), alpha = 0.5, outlier.shape = NA) +
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

### Institution
data_prop_LFS <- data_prop_LFS %>%
  group_by(institution) %>%
  arrange(p_20_150)
fig_inst <- ggplot(data_prop_LFS) +
  geom_boxplot(aes(institution, p_20_150), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(institution, p_20_150, group = institution, color = cancer), alpha = 0.75,
             position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_inst, aes(x = institution, y = 0.03, label = N), size = 4) +
  geom_hline(yintercept = cohort_median, linetype = "dashed", color = "red", size = 0.25) +
  ggtitle("Institution") + 
  xlab("") +
  ylab("Proportion under 150bp") +
  scale_fill_manual(values = c("grey65", "#FB9A99")) +
  scale_color_manual(labels = c("Negative", "Positive"),
                     values = c("black", "red")) +
  labs(fill = "Patient", color = "Cancer Status") +
  theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_y_continuous(limits=c(0, 0.45), expand = c(0,0))
fig_inst

### plot read groups
data_reads <- data_reads %>%
  group_by(V1) %>%
  arrange(p_20_150)
fig_reads <- ggplot(data_reads) +
  geom_boxplot(aes(V1, p_20_150), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(V1, p_20_150, group = V1, color = cancer), alpha = 0.75,
             position = position_dodge2(.5), size = 1.5, pch = 16) +
  geom_text(data = data_stats_group, aes(x = V1, y = 0.03, label = N), size = 4) +
  geom_text(data = data_stats_group, aes(x = V1, y = 0.41, label = annot), size = 5) +
  geom_hline(yintercept = cohort_median, linetype = "dashed", color = "red", size = 0.25) +
  ggtitle("Read Group") + 
  xlab("") +
  ylab("Proportion under 150bp") +
  scale_fill_manual(values = c("grey65", "#FB9A99")) +
  scale_color_manual(labels = c("Healthy", "Negative", "Positive"),
                     values = c("blue", "black", "red")) +
  labs(fill = "Patient", color = "Cancer Status") +
  theme +
  theme(legend.position = "right",
        axis.text.x = element_blank()) +
  scale_y_continuous(limits=c(0, 0.45), expand = c(0,0))
fig_reads

fig_cluster <- ggplot(clusters) +
  geom_bar(aes(run.ID, X), stat = "identity") +
  ggtitle("") + 
  xlab("Sequencing Run") +
  ylab("Median Flow Cell\nCluster Density") +
  theme
fig_cluster

fig_read_groups <- fig_reads/fig_cluster
fig_read_groups <- fig_read_groups + plot_layout(heights = c(3, 1))
fig_read_groups

### Plot regression of cluster densities
fig_reg <- ggplot(data_reads) +
  geom_point(aes(X, p_20_150, color = cancer), alpha = 0.5) +
  geom_smooth(aes(X, p_20_150, color = cancer), method = "lm", se = FALSE, size = 1, alpha = 0.5) +
  stat_regline_equation(label.y = c(0, 0.05, 0.1), label.x = 100, aes(X, p_20_150, color = cancer, label = ..rr.label..), size = 5) +
  ggtitle("Regression") + 
  xlab("Median Cluster Density") +
  ylab("Proportion under 150bp") +
  theme +
  scale_color_manual(labels = c("Healthy", "Negative", "Positive"),
                     values = c("blue", "black", "red")) +
  labs(fill = "Patient", color = "Cancer Status")
fig_reg

ggsave(file.path(outdir, "fragment_proportions_lfs.pdf"), fig, device = "pdf", width = 2, height = 4, units = "in")
ggsave(file.path(outdir, "fragment_proportions_lfsh.pdf"), fig_neg, device = "pdf", width = 2, height = 4, units = "in")
ggsave(file.path(outdir, "fragment_proportions_lfs_status.pdf"), fig_lfs, device = "pdf", width = 2.5, height = 4, units = "in")
ggsave(file.path(outdir, "fragment_proportions_lfs_history.pdf"), fig_previous, device = "pdf", width = 2, height = 4, units = "in")
ggsave(file.path(outdir, "fragment_proportions_age_reg.pdf"), fig_age, device = "pdf", width = 4, height = 4, units = "in")
ggsave(file.path(outdir, "fragment_proportions_type.pdf"), fig_type, device = "pdf", width = 3.5, height = 4.5, units = "in")
ggsave(file.path(outdir, "fragment_proportions_germ.pdf"), fig_germ, device = "pdf", width = 10, height = 5, units = "in")
ggsave(file.path(outdir, "fragment_proportions_fam.pdf"), fig_fam, device = "pdf", width = 8, height = 5, units = "in")
ggsave(file.path(outdir, "fragment_proportions_institution.pdf"), fig_inst, device = "pdf", width = 2, height = 4, units = "in")
ggsave(file.path(outdir, "fragment_proportions_read_group.pdf"), fig_read_groups, device = "pdf", width = 10, height = 6, units = "in")
ggsave(file.path(outdir, "fragment_proportions_read_reg.pdf"), fig_reg, device = "pdf", width = 6, height = 6, units = "in")


