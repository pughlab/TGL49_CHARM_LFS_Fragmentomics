library(dplyr)
library(matrixStats)
library(psych)
library(ComplexHeatmap)

### Set variables
path <- "data/fragment_ratio"
outdir <- ""
healthy <- "hbc/fragment_ratio"


### Import data (Starting with the 5Mb ratios)
de_table <- read.delim(file.path(path, "fragment_coverage_differential.txt"))
data_normal <- read.delim(file.path(healthy, "TGL49_HBC_frags_100kb.txt"))
data_cov <- read.delim(file.path(path, "CHARM_LFS_frags_100kb.txt"))
samples <- read.delim("sample_list.txt")

### Format healthy controls and calculate median
row.names(data_normal) <- with(data_normal, paste0(seqnames, "_", start))
data_normal <- data_normal[, -c(1:3)]
sums <- colSums(data_normal, na.rm = TRUE)
data_normal <- sweep(data_normal, 2, sums, "/")
data_normal <- data_normal*1000000

median_hbc <- rowMedians(as.matrix(data_normal), na.rm = TRUE)

### Calculate medians for each cluster
row.names(data_cov) <- with(data_cov, paste0(seqnames, "_", start))
data_cov <- data_cov[, -c(1:3)]
sums <- colSums(data_cov, na.rm = TRUE)
data_cov <- sweep(data_cov, 2, sums, "/")
data_cov <- data_cov*1000000

median_cohort <- rowMedians(as.matrix(data_cov), na.rm = TRUE)
median_lof <- rowMedians(as.matrix(data_cov[,colnames(data_cov) %in% samples$sWGS[samples$mutation_type == "LOF"]]), na.rm = TRUE)
median_splice <- rowMedians(as.matrix(data_cov[,colnames(data_cov) %in% samples$sWGS[samples$mutation_type == "Splice"]]), na.rm = TRUE)
median_c2 <- rowMedians(as.matrix(data_cov[,colnames(data_cov) %in% samples$sWGS[samples$mutation_type == "2"]]), na.rm = TRUE)
median_c3 <- rowMedians(as.matrix(data_cov[,colnames(data_cov) %in% samples$sWGS[samples$mutation_type == "3"]]), na.rm = TRUE)

### Make correlation table
de_corr <- data.frame(median_hbc, median_cohort, median_c2, median_c3, median_lof, median_splice)
row.names(de_corr) <- row.names(data_cov)
colnames(de_corr) <- c("HBC", "cohort", "c2", "c3", "lof", "splice")

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
de_sum_melt$variable <- factor(de_sum_melt$variable, levels = c("Cohort", "Missense 3", "Splice", "Missense 2", "LOF"))
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
ann_text$type <- "LOF"
ann_text$type <- factor(ann_text$type, levels = c("Cohort", "Missense 3", "Splice", "Missense 2", "LOF"))

### Plot correlation table
pdf(file.path(outdir, "fragment_coverage_correlation_table.pdf"), height = 3, width = 3)
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
ggsave(file.path(outdir, "fragment_coverage_chromosomes.pdf"), fig_chr, device = "pdf", width = 10, height = 4.5, units = "in")

