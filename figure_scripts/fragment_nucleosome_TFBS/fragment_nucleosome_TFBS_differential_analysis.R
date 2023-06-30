library(dplyr)
library(matrixStats)
library(lemon)
library(scales)
library(ggrepel)

### Set paths
outdir <- ""

data <- read.delim("data/griffin/fragment_nucleosome_TFBS_differential.txt")
data_samples <- read.delim("sample_list.txt")

### Seperate out features
data$site <- gsub("_.*","", data$feature)
data$variable <- gsub(".*_","", data$feature)

### Format plotting table
data$variable <- factor(data$variable, levels = c("midpoint", "amp"),
                        labels = c("Midpoint\nCoverage", "Amplitude"))

data_change <- data[, colnames(data) %in% c("variable", "site", "cohort", "negative", "positive", "cancer_status")]
data_change <- reshape2::melt(data_change, id = c("variable", "site"))
colnames(data_change) <- c("feature", "site", "type", "value")
data_p <- data[, colnames(data) %in% c("variable", "site", "padj_cohort", "padj_negative", "padj_positive", "padj_cancer_status")]
data_p <- reshape2::melt(data_p, id = c("variable", "site"))
colnames(data_p) <- c("feature", "site", "type", "value")
data_p$type <- gsub("padj_", "", data_p$type)

data <- merge(data_change, data_p, by = c("feature", "site", "type"))
colnames(data) <- c("feature", "site", "type", "change", "padj")

data$type <- factor(data$type, levels = c("cohort", "negative", "positive", "cancer_status"),
                    labels = c("LFS vs Healthy", "LFS Cancer Negative\nvs Healthy", "LFS Cancer Positive\nvs Healthy", "LFS Cancer\nPositive vs Negative"))

### Set color annotation
data$color <- "none"
data$color <- ifelse(data$change > 1 & data$padj < 0.01, "positive", data$color)
data$color <- ifelse(data$change < -1 & data$padj < 0.01, "negative", data$color)

### Calculate stats
stats <- data %>%
  group_by(feature, type, color) %>%
  dplyr::summarise(N=n())

### Split tables
data_cohort <- data[data$type == "LFS vs Healthy", ]
data <- data[!(data$type == "LFS vs Healthy"), ]

### Make annotation table
data_cohort$annotation <- ifelse(data_cohort$feature == "Midpoint\nCoverage" & data_cohort$site == "TP53", data_cohort$site, "")

### Set y-axis function
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

### Plot Volcano plots
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
               axis.title = element_text(size = 12),
               axis.line = element_line(colour = "black"),
               axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
               axis.text = element_text(size = 10),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.background = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.key = element_rect(fill = "white"),
               legend.spacing.y = unit(0, "mm"),
               legend.key.size = unit(3, "mm"),
               strip.background = element_rect(color = NA, fill = NA),
               strip.text = element_text(face = "bold", size = 10))

plot <- ggplot(data, aes(change, padj, color = color)) +
  geom_point(size = 0.75, pch = 16) +
  geom_hline(yintercept = 0.01, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) +
  facet_grid(type ~ feature, scales = "free") +
  facet_rep_grid(type ~ feature, scales = "free") +
  xlab("Median Distance from Comparator (SDs)") + 
  ylab("p-adjusted") +
  scale_color_manual(values = c("none" = "black", "positive" = "red", "negative" = "blue")) +
  theme +
  scale_y_continuous(trans=reverselog_trans(10))
plot

ggsave(file.path(outdir, paste0("fragment_nucleosome_TFBS_volcano_all.pdf")), plot, device = "pdf", width = 6, height = 7, units = "in")

plot2 <- ggplot(data_cohort, aes(change, padj)) +
  geom_point(aes(color = color), size = 0.75, pch = 16) +
  geom_hline(yintercept = 0.01, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) +
  geom_text_repel(aes(label = annotation), size = 3, max.overlaps = Inf, xlim = c(-1,0), ylim = c(6,7)) +
  facet_wrap(.~feature, scales = "free", ncol = 2) +
  xlab("Median Distance (SDs)") + 
  ylab("p-adjusted\n(LFS vs Healthy)") +
  scale_color_manual(values = c("none" = "black", "positive" = "red", "negative" = "blue")) +
  theme +
  scale_y_continuous(trans=reverselog_trans(10))
plot2

ggsave(file.path(outdir, paste0("fragment_nucleosome_TFBS_volcano_cohort.pdf")), plot2, device = "pdf", width = 6, height = 3, units = "in")

### Find TFs that are DE
data_de <- data[data$feature == "Midpoint\nCoverage" &
                  !(data$type == "LFS Cancer\nPositive vs Negative"), ]

df <- data_de %>%
  group_by(feature, site) %>%
  dplyr::mutate(Diff = change - lag(change))
df <- df[complete.cases(df), ]

data_de <- merge(data_de, df[, c("site", "Diff")], by = "site")

data_de$direction <- ifelse(data_de$Diff > 0.5, "Decreased\nAccessibility", 
                            ifelse(data_de$Diff < -0.5, "Increased\nAccessibility", "No Change"))
data_de <- data_de[data_de$direction %in% c("Decreased\nAccessibility", "Increased\nAccessibility"), ]
data_de$type <- factor(data_de$type, levels = c("LFS Cancer Negative\nvs Healthy", "LFS Cancer Positive\nvs Healthy"),
                       labels = c("Negative", "Positive"))
data_de$annotation <- ifelse((data_de$Diff > 0.5 | data_de$Diff < -0.5) &
                               data_de$type == "Positive", data_de$site, "")
data_de$color <- ifelse(data_de$direction == "Decreased\nAccessibility", "red",
                        ifelse(data_de$direction == "Increased\nAccessibility", "blue", "none"))

data_de$annotation <- ifelse(data_de$annotation == "", NA, data_de$annotation)

### Format Annotations
data_de$annotation <- ifelse(data_de$annotation %in% c("ETV2", "GRHL2", "GATA3", "HES1", "POU5F1"),
                             data_de$annotation, NA)
data_de$annotation <- ifelse(data_de$type == "Positive" & 
                               data_de$direction == "No Change" &
                               (data_de$change > 1.5 | data_de$change < -1.5), 
                             data_de$site, data_de$annotation)

labels <- data_de[complete.cases(data_de), ]

### Plot changes
plot_de <- ggplot(data_de, aes(type, change)) +
  geom_point(aes(color = color), size = 3, pch = 16, alpha = 0.5) +
  geom_line(aes(group = site), alpha = 0.5) +
  geom_text_repel(data = labels, aes(label = annotation), size = 2.5, xlim = c(2.5,3), force = 20, max.overlaps = 20) +
  facet_grid(.~direction, scales = "free") +
  facet_rep_grid(.~direction, scales = "free") +
  xlab("Cancer Status") + 
  ylab("Median Distance vs Healthy (SDs)") +
  scale_color_manual(values = c("blue", "red")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = "white"),
        legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(3, "mm"),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(face = "bold", size = 10)) +
  scale_x_discrete(expand = expansion(mult=c(0.5,1.5)))
plot_de

ggsave(file.path(outdir, paste0("fragment_nucleosome_TFBS_cancer_status_diff.pdf")), plot_de, device = "pdf", width = 4, height = 3, units = "in")

de_list <- data_de[data_de$direction %in% c("Increased\nAccessibility", "Decreased\nAccessibility"), ]
de_list$direction <- ifelse(de_list$direction == "Increased\nAccessibility", "Increase", "Decrease")
de_list <- de_list[!(duplicated(de_list$site)), ]
de_list <- de_list[, c("site", "direction")]
de_list <- de_list[order(de_list$direction), ]
write.table(de_list, file.path(outdir, "fragment_nucleosome_TFBS_delist.txt"), sep = "\t", row.names = FALSE)

stats_lfs <- as.data.table(table(de_list$direction))
