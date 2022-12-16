library(tidyverse)
library(dplyr)
library(matrixStats)

### Set paths
path <- "/Volumes/GoogleDrive/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/dnase1l3"

### Import data 
cna <- read.delim(list.files(path, "cna", full.names = TRUE))
mut <- read.delim(list.files(path, "mut", full.names = TRUE))
svs <- read.delim(list.files(path, "structural", full.names = TRUE))
rna <- read.delim(list.files(path, "TCGA_expression", full.names = TRUE))

### Format data
rna <- rna[complete.cases(rna), ]
rna <- rna[!(rna$sample_type %like% "Metastatic" |
               rna$sample_type %like% "Additional" |
               rna$sample_type %like% "Recurrent"), ]
rna$DNASE1L3 <- 2^rna$DNASE1L3 - 0.001

mut$TP53 <- ifelse(mut$TP53 == "NS", NA, mut$TP53)
mut <- mut[complete.cases(mut), ]
TP53 <- tidyr::separate(data = mut, col = TP53, sep = " ", into = c("TP53_1", "TP53_2", "TP53_3"), remove = FALSE)
TP53[is.na(TP53)] <- 0
TP53$TP53_1 <- ifelse(TP53$TP53_1 == "WT", 0, 1)
TP53$TP53_2 <- ifelse(TP53$TP53_2 == 0, 0, 1)
TP53$TP53_3 <- ifelse(TP53$TP53_3 == 0, 0, 1)
TP53$sum <- rowSums(TP53[,7:9])

cna$TP53 <- ifelse(cna$TP53 == "NP", NA, cna$TP53)
cna <- cna[complete.cases(cna), ]
cna$TP53 <- ifelse(cna$TP53 == 0, 0, 1)

svs[is.na(svs)] <- 0
svs$TP53 <- ifelse(svs$TP53 == 0, 0, 1)

### Merge data together
data <- merge(cna[, c(1,2,6)], TP53[, c(2,10)], by = "SAMPLE_ID")
data <- merge(data, svs[, c(2,6)], by = "SAMPLE_ID")
data$variants <- rowSums(data[, 3:5])
data$status <- ifelse(data$variants == 0, "Wildtype",
                      ifelse(data$variants == 1, "Heterozygous", "Homozygous"))
data <- merge(data[, c(1,2,7)], rna, by.x = "SAMPLE_ID", by.y = "sample", all = TRUE)
data <- data[!(is.na(data$DNASE1L3)), ]
data$status <- ifelse(data$sample_type %like% "Normal", "Normal", data$status)
data <- data[!(is.na(data$status) |
                 data$cancer.type.abbreviation == ""), ]

### Format data
names(data)[names(data) == "cancer.type.abbreviation"] <- "type"
data$status <- factor(data$status, levels = c("Wildtype", "Heterozygous", "Homozygous", "Normal"))

### Calculate stats
data_stats <- data %>%
  group_by(status, type) %>%
  dplyr::summarise(mean = mean(DNASE1L3),
                   median = median(DNASE1L3),
                   sd = sd(DNASE1L3),
                   N=n())
data_stats <- data_stats[data_stats$N > 1, ]

remove <- c("PCPG", "TGCT", "THCA", "UVM", "SKCM")
data_stats <- data_stats[!(data_stats$type %in% remove), ]
data <- data[!(data$type %in% remove), ]

types <- unique(data_stats$type)
status <- unique(data_stats$status)
t_test <- c()

for (stat in status) {
  for (type in types) {
    if(stat == "Normal" & type %in% c("ACC", "DLBC", "GBM", "LGG", "LAML", "MESO", "OV", "UCS")) next
    a <- t.test(data$DNASE1L3[data$type == type & data$status == "Wildtype"],
                data$DNASE1L3[data$type == type & data$status == stat])$p.value
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
               legend.position = "bottom",
               legend.key = element_rect(fill = "white"),
               legend.title = element_text(size = 12),
               legend.text = element_text(size = 12),
               strip.background = element_blank(),
               strip.text = element_text(size = 13),
               axis.text = element_text(size = 13),
               axis.text.x = element_blank(), 
               axis.title = element_text(size = 13))

### Plot curves
fig <- ggplot(data) +
  geom_boxplot(aes(status, DNASE1L3, fill = status)) +
  geom_text(data = data_stats, aes(status, y = 400, label = annot), size = 4) +
  geom_text(data = data_stats, aes(status, y = 0.005, label = N), size = 2.5) +
  xlab("Cohort") + 
  ylab("mRNA Expression (RSEM)") +
  labs(color = "", fill = "") +
  facet_wrap(.~type, nrow = 2) +
  scale_fill_manual(values = c("grey65", "#67a9cf", "#ef8a62", "#B2DF8A")) +
  ggtitle("DNASE1L3 Expression") + 
  scale_y_log10() +
  #scale_y_continuous(limits = c(0, 100)) +
  theme
fig

ggsave(file.path(path, "DNASE1L3_expression.pdf"), fig, width = 12, height = 5)






