library(dplyr)
library(matrixStats)
library(psych)
library(ComplexHeatmap)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_coverage"
healthy <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC"

### Import data (Starting with the 5Mb ratios)
data_cov <- read.delim(list.files(file.path(path, "fragmentomics"), "frag", full.names = TRUE))
samples <- read.delim(file.path(path, "samples/sample_list.txt"))
data_normal <- read.delim(list.files(file.path(healthy, "fragmentomics"), "frag", full.names = TRUE))

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]

### Keep only plasma samples and order based on clinical information
samples <- samples[!(samples$germline_mutation == "unknown"), ]
samples <- samples[samples$sWGS %in% colnames(data_cov), ]
samples <- samples[samples$cancer_status == "negative", ]

### Set chromosomes
data_chr <- data_cov[, c("seqnames", "start", "end")]
data_chr$name <- with(data_chr, paste0(seqnames, "_", start))

### Calculate the healthy median
row.names(data_normal) <- with(data_normal, paste0(seqnames, "_", start))
data_normal <- as.matrix(data_normal[, -c(1:3)])
sums <- colSums(data_normal, na.rm = TRUE)
data_normal <- sweep(data_normal, 2, sums, "/")
data_normal <- data_normal*1000000
data_normal <- data_normal[complete.cases(data_normal), ]

normal_median <- rowMedians(data_normal, na.rm = TRUE)
normal_sd <- rowSds(data_normal, na.rm = TRUE)

### Make matrix for healthy median
normal_medians <- as.matrix(normal_median)
row.names(normal_medians) <- row.names(data_normal)
normal_sds <- as.matrix(normal_sd)
row.names(normal_sds) <- row.names(data_normal)

### Format ratio sheet and order samples properly
row.names(data_cov) <- with(data_cov, paste0(seqnames, "_", start))
data_cov <- data_cov[, samples$sWGS]
sums <- colSums(data_cov, na.rm = TRUE)
data_cov <- sweep(data_cov, 2, sums, "/")
data_cov <- data_cov*1000000
data_cov <- data_cov[row.names(data_normal), ]

### Perform differential analysis (Format samples and data)
data_de <- bind_cols(data_cov, data_normal)
data_de <- as.data.frame(t(data_de))

samples <- samples[samples$cancer_status == "negative", ]
samples_de <- as.data.frame(cbind(samples[, c("sWGS", "mutation_type")], "LFS"))
colnames(samples_de) <- c("sWGS", "type", "diag")
normal_de <- as.data.frame(cbind(colnames(data_normal), "HBC", "HBC"))
colnames(normal_de) <- c("sWGS", "type", "diag")
samples_de <- bind_rows(samples_de, normal_de)

data_de <- merge(data_de, samples_de, by.x = "row.names", by.y = "sWGS")
data_de_melt <- reshape2::melt(data_de, id = c("Row.names", "diag", "type"))

### Calculate t-tests
t_test <- c()
vars <- unique(data_de_melt$variable)
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$diag == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$diag == "LFS" & data_de_melt$variable == var])$p.value
  t_test <- c(t_test, a)
}

t_test2 <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "2" & data_de_melt$variable == var])$p.value
  t_test2 <- c(t_test2, a)
}

t_test3 <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "3" & data_de_melt$variable == var])$p.value
  t_test3 <- c(t_test3, a)
}

t_testlof <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "LOF" & data_de_melt$variable == var])$p.value
  t_testlof <- c(t_testlof, a)
}

t_testsp <- c()
for (var in vars) {
  a <- t.test(data_de_melt$value[data_de_melt$type == "HBC" & data_de_melt$variable == var], 
              data_de_melt$value[data_de_melt$type == "Splice" & data_de_melt$variable == var])$p.value
  t_testsp <- c(t_testsp, a)
}

t_test <- as.data.frame(cbind(as.character(vars), t_test, t_test2, t_test3, t_testlof, t_testsp))
t_test$t_test <- as.numeric(t_test$t_test)
t_test$t_test2 <- as.numeric(t_test$t_test2)
t_test$t_test3 <- as.numeric(t_test$t_test3)
t_test$t_testlof <- as.numeric(t_test$t_testlof)
t_test$t_testsp <- as.numeric(t_test$t_testsp)

### Calculate the distance from the healthy median
data_change <- (data_cov - normal_median)/normal_sd
data_normal <- (data_normal - normal_median)/normal_sd

### Set Median change of LFS cohort from HBCs and within mutation groups
data_freq <- rowMedians(as.matrix(data_change), na.rm = TRUE)
data_freq <- as.matrix(data_freq)
row.names(data_freq) <- row.names(data_change)

data_c2 <- data_change[, colnames(data_change) %in% samples$sWGS[samples$mutation_type == 2]]
data_freq2 <- rowMedians(as.matrix(data_c2), na.rm = TRUE)
data_freq2 <- as.matrix(data_freq2)
row.names(data_freq2) <- row.names(data_change)

data_c3 <- data_change[, colnames(data_change) %in% samples$sWGS[samples$mutation_type == 3]]
data_freq3 <- rowMedians(as.matrix(data_c3), na.rm = TRUE)
data_freq3 <- as.matrix(data_freq3)
row.names(data_freq3) <- row.names(data_change)

data_lof <- data_change[, colnames(data_change) %in% samples$sWGS[samples$mutation_type == "LOF"]]
data_freqlof <- rowMedians(as.matrix(data_lof), na.rm = TRUE)
data_freqlof <- as.matrix(data_freqlof)
row.names(data_freqlof) <- row.names(data_change)

data_sp <- data_change[, colnames(data_change) %in% samples$sWGS[samples$mutation_type == "Splice"]]
data_freqsp <- rowMedians(as.matrix(data_sp), na.rm = TRUE)
data_freqsp <- as.matrix(data_freqsp)
row.names(data_freqsp) <- row.names(data_change)

data_freq <- bind_cols(data_freq, data_freq2, data_freq3, data_freqlof, data_freqsp)
data_freq <- as.data.frame(data_freq)
row.names(data_freq) <- row.names(data_change)
colnames(data_freq) <- c("all", "c2", "c3", "lof", "splice")

### Make Differential table and plot
de_table <- merge(data_freq, t_test, by.x = "row.names", by.y = "V1")
de_table$chr <- gsub("_.*", "", de_table$Row.names)
de_table <- de_table[order(factor(de_table$Row.names, levels = data_chr$name)), ]
de_table <- bind_cols(de_table, normal_medians)
de_table <- bind_cols(de_table, normal_sds)
colnames(de_table) <- c("location", "cohort", "c2", "c3", "lof", "splice", "padj_cohort", "padj_c2", 
                        "padj_c3", "padj_lof", "padj_splice", "chr", "healthy_median", "healthy_sd")
write.table(de_table, file.path(outdir, "fragment_coverage_differential.txt"), sep = "\t", row.names = FALSE)
