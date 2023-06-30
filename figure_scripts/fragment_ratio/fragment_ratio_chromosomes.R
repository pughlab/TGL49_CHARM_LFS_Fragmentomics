library(dplyr)
library(matrixStats)
library(psych)
library(ComplexHeatmap)
library(circlize)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_ratio"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"

### Import data (Starting with the 5Mb ratios)
de_table <- read.delim(file.path(outdir, "fragment_ratio_differential.txt"))
data_normal <- read.delim(file.path(healthy, "fragmentomics", "TGL49_HBC_ratio_5Mb.txt"))

### Calculate # and stats
de_up <- de_table[de_table$padj_cohort < 0.05 & de_table$cohort > 0, ]
de_up_median <- median((abs(de_up$healthy_median) + (de_up$healthy_sd*de_up$cohort))/abs(de_up$healthy_median))

de_down <- de_table[de_table$padj_cohort < 0.05 & de_table$cohort < 0, ]
de_down_median <- median((abs(de_down$healthy_median) + (de_down$healthy_sd*de_down$cohort))/abs(de_down$healthy_median))

### Format healthy controls and calculate a mean distance from median
row.names(data_normal) <- with(data_normal, paste0(seqnames, "_", start))
data_normal <- as.matrix(data_normal[, -c(1:4)])

data_normal <- (data_normal - de_table$healthy_median)
HBC <- as.matrix(rowMeans2(data_normal))
row.names(HBC) <- row.names(data_normal)

### Make correlation table
de_corr <- de_table[, c("cohort", "c2", "c3", "lof", "splice")]
row.names(de_corr) <- de_table$location
de_corr <- bind_cols(de_corr, HBC)
names(de_corr)[names(de_corr) == "...6"] <- "HBC"
de_corr <- de_corr[row.names(data_normal) , c("HBC", "cohort", "c2", "c3", "lof", "splice")]

types <- colnames(de_corr)
y <- c("type", types)
for (type in types) {
  x <- c(type)
  for (type2 in types) {
    a <- corr.test(de_corr[, type], de_corr[, type2])[["r"]]
    x <- c(x, a)
  }
  y <- rbind(y, x)
}

### Format correlation table and set graphics
corr_table <- as.data.frame(y)
corr_table <- corr_table[-1, -1]
colnames(corr_table) <- c("Healthy", "Cohort", "Missense 2", "Missense 3", "LOF", "Splice")
row.names(corr_table) <- c("Healthy", "Cohort", "Missense 2", "Missense 3", "LOF", "Splice")
corr_table <- mutate_all(corr_table, function(x) as.numeric(as.character(x)))
corr_table <- as.matrix(corr_table)

col_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), 
                      c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", corr_table[i, j]), x, y, gp = gpar(fontsize = 8))
}

### Convert to binaries
de_table$cohort <- ifelse(de_table$padj_cohort < 0.05 & de_table$cohort < 0, -1,
                          ifelse(de_table$padj_cohort < 0.05 & de_table$cohort > 0, 1, 0))
de_table$c2 <- ifelse(de_table$padj_c2 < 0.05 & de_table$c2 < 0, -1,
                       ifelse(de_table$padj_c2 < 0.05 & de_table$c3 > 0, 1, 0))
de_table$c3 <- ifelse(de_table$padj_c3 < 0.05 & de_table$c3 < 0, -1,
                       ifelse(de_table$padj_c3 < 0.05 & de_table$c3 > 0, 1, 0))
de_table$lof <- ifelse(de_table$padj_lof < 0.05 & de_table$lof < 0, -1,
                       ifelse(de_table$padj_lof < 0.05 & de_table$lof > 0, 1, 0))
de_table$splice <- ifelse(de_table$padj_splice < 0.05 & de_table$splice < 0, -1,
                       ifelse(de_table$padj_splice < 0.05 & de_table$splice > 0, 1, 0))

de_sum <- de_table[, c("chr", "cohort", "c2", "c3", "lof", "splice")]
de_sum <- de_sum %>%
  group_by(chr) %>%
  dplyr::summarise(all_pos = sum(cohort[cohort == 1]),
            all_neg = sum(cohort[cohort == -1]),
            c2_pos = sum(c2[c2 == 1]),
            c2_neg = sum(c2[c2 == -1]),
            c3_pos = sum(c3[c3 == 1]),
            c3_neg = sum(c3[c3 == -1]),
            lof_pos = sum(lof[lof == 1]),
            lof_neg = sum(lof[lof == -1]),
            sp_pos = sum(splice[splice == 1]),
            sp_neg = sum(splice[splice == -1]),
            count = n())
cohort_neg <- sum(de_sum$all_neg)
cohort_pos <- sum(de_sum$all_pos)
chromosomes <- paste0("chr", c(1:22))
de_sum <- de_sum[order(factor(de_sum$chr, levels = chromosomes)), ]
count <- de_sum$count
de_sum <- de_sum[ , c(2:11)]/de_sum$count
de_sum$chr <- chromosomes
de_sum$count <- count
de_sum$chr <- factor(de_sum$chr, levels = chromosomes)

### Melt and format summed data for plotting
de_sum_melt <- reshape2::melt(de_sum, id = c("chr", "count"))
de_sum_melt$variable <- factor(de_sum_melt$variable, levels = c("all_pos", "all_neg", "c2_pos", "c2_neg", "c3_neg", "c3_pos", "lof_pos", "lof_neg", "sp_pos", "sp_neg"),
                               labels = c("Cohort", "Cohort", "Missense 2", "Missense 2", "Missense 3", "Missense 3", "LOF", "LOF", "Splice", "Splice"))
de_sum_melt$variable <- factor(de_sum_melt$variable, levels = c("Cohort", "Splice", "LOF", "Missense 3", "Missense 2"))
de_sum_melt_pos <- de_sum_melt[de_sum_melt$value >= 0, ]
de_sum_melt_neg <- de_sum_melt[de_sum_melt$value <= 0, ]
de_sum_melt <- merge(de_sum_melt_pos, de_sum_melt_neg, by = c("chr", "variable", "count"), all = TRUE)
de_sum_melt[is.na(de_sum_melt)] <- 0
colnames(de_sum_melt) <- c("chr", "type", "count", "pos", "neg")
de_sum_melt$chr <- factor(de_sum_melt$chr, levels = chromosomes)

### Make text annotation
ann_text <- de_sum_melt %>%
  group_by(chr) %>%
  dplyr::summarise(count = unique(count))
ann_text$type <- "Missense 2"
ann_text$type <- factor(ann_text$type, levels = c("Cohort", "Missense 3", "Splice", "Missense 2", "LOF"))

### Perform Chi-squared tests
data_chi <- de_sum_melt
data_chi$pos <- data_chi$count*data_chi$pos
data_chi$neg <- -data_chi$count*data_chi$neg
data_chi$none <- data_chi$count - data_chi$pos - data_chi$neg
data_chi <- data_chi[data_chi$type == "Cohort", ]
data_chi <- data_chi[!(data_chi$pos == 0 & data_chi$neg == 0), ]

stats_chi_chr <- data.frame()
for (i in data_chi$chr) {
    x <- data_chi[data_chi$chr == i, c("pos", "neg", "none")]
    x <- colSums(x)
    y <- data_chi[!(data_chi$chr == i), c("pos", "neg", "none")]
    y <- colSums(y)
    z <- data.frame(a = x,
                    b = y)
    p <- c(i, chisq.test(z)$p.value)
    stats_chi_chr <- rbind(stats_chi_chr, p)
}

data_chi <- de_sum_melt
data_chi$pos <- data_chi$count*data_chi$pos
data_chi$neg <- -data_chi$count*data_chi$neg
data_chi$none <- data_chi$count - data_chi$pos - data_chi$neg
data_chi <- data_chi %>%
  group_by(type) %>%
  dplyr::summarise(up = sum(pos),
                   down = sum(neg))
data_chi$none <- 512 - data_chi$up - data_chi$down
names <- data_chi$type
data_chi <- as.data.frame(t(data_chi[, 2:4]))
colnames(data_chi) <- names

stats_chi_class <- c()
for (i in names) {
  x <- data_chi[, colnames(data_chi) == i]
  y <- data_chi[, colnames(data_chi) == "Cohort"]
  z <- cbind(x, y)
  p <- chisq.test(z)$p.value
  stats_chi_class <- c(stats_chi_class, p)
}
stats_chi_class <- data.frame(class = names,
                              pvalue = stats_chi_class)

### Plot correlation table
pdf(file.path(outdir, "fragment_genome_correlation_table.pdf"), height = 3, width = 3)
fig_corr <- Heatmap(corr_table,
                    col = col_fun,
                    cell_fun = cell_fun,
                    show_heatmap_legend = FALSE,
                    row_names_gp = gpar(fontsize = 10),
                    column_names_gp = gpar(fontsize = 10),
                    height = unit(0.75*nrow(corr_table), "cm"),
                    width = unit(0.75*nrow(corr_table), "cm"))
draw(fig_corr)
dev.off()

### Plot summed data
fig_chr <- ggplot(de_sum_melt, aes(x = forcats::fct_rev(chr))) +
  geom_bar(aes(y = pos, fill = "red"), stat = "identity") +
  geom_bar(aes(y = neg, fill = "blue"), stat = "identity") +
  geom_text(data = ann_text, aes(x = chr, y = 0.8, label = count), size = 4) +
  facet_wrap(.~type, nrow = 1) +
  scale_fill_manual(values = c("red" = "#fb9a99", "blue" = "#a6cee3")) +
  labs(fill = "") +
  xlab("Chromosome") + 
  ylab("Fraction of bins") +
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5, size = 13), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0)) +
  scale_y_continuous(limits=c(-1, 1), expand = c(0,0), breaks = c(-0.5, 0, 0.5)) + 
  coord_flip()
fig_chr  
ggsave(file.path(outdir, "fragment_genome_chromosomes.pdf"), fig_chr, device = "pdf", width = 10, height = 4.5, units = "in")

