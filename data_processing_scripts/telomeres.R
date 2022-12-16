library(tidyverse)
library(plyr)
library(matrixStats)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/telomerehunter/output"
path2 <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/telseq/output"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/telomeres"
project <- "CHARM_LFS"

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in data
samples <- read.delim(samples)

filenames <- list.files(path = path, pattern = "_TVRs.txt", recursive = TRUE, full.names = TRUE)
names <- list.files(path = path, pattern = "_TVRs.txt", recursive = TRUE, full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

datalist <- lapply(filenames, function(x){read.delim(file = x)[,1:2]})
contexts <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Pattern", all = TRUE), datalist)
contexts <- as.data.frame(contexts)

filenames <- list.files(path = path, pattern = "_tumor_gc_content.tsv", recursive = TRUE, full.names = TRUE)
filenames <- filenames[!(grepl("intratelomeric", filenames))]
datalist <- lapply(filenames, function(x){read.delim(file = x)})
reads_gc <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gc_content_percent", all = TRUE), datalist)
reads_gc <- as.data.frame(reads_gc)

filenames <- list.files(path = path, pattern = "_summary.tsv", recursive = TRUE, full.names = TRUE)
filenames <- filenames[!(grepl("tumor", filenames))]
datalist <- lapply(filenames, function(x){read.delim(file = x)})
summaries <- do.call(rbind, datalist)
summaries <- as.data.frame(summaries)

filenames <- list.files(path = path2, pattern = "telseq.txt", recursive = TRUE, full.names = TRUE)
datalist <- lapply(filenames, function(x){read.delim(file = x, skip = 2)})
telseq <- do.call(rbind, datalist)
telseq <- as.data.frame(telseq)

### Format TelomereHunter outputs
# telomere contexts
row.names(contexts) <- contexts$Pattern
contexts <- contexts[, -1]
colnames(contexts) <- names
contexts[is.na(contexts)] <- 0
contexts <- t(contexts)
contexts <- as.data.frame(contexts)
contexts$sums <- rowSums2(as.matrix(contexts))

# GC Reads
row.names(reads_gc) <- reads_gc$gc_content_percent
reads_gc <- reads_gc[, -1]
colnames(reads_gc) <- names
reads_gc <- reads_gc[row.names(reads_gc) %in% c(48, 49, 50, 51, 52), ]
reads_gc <- t(reads_gc)
reads_gc <- as.data.frame(reads_gc)
reads_gc$total <- rowSums2(as.matrix(reads_gc))

# Summaries
summaries <- summaries[, colnames(summaries) %in% c("PID", "total_reads", "tel_reads", "intratel_reads", "total_reads_with_tel_gc", "tel_content")]
colnames(summaries) <- c("sample", "total_reads", "tel_reads", "intraltel_reads", "reads_gc", "tel_content")

### Calculate normalized fragments and percent
contexts_norm <- (contexts*1000000)/reads_gc$total
contexts_per <- contexts_norm/contexts_norm$sums

contexts_norm$sample <- row.names(contexts_norm)
contexts_norm <- contexts_norm %>%
  dplyr::select(sample, everything())

contexts_per$sample <- row.names(contexts_per)
contexts_per <- contexts_per %>%
  dplyr::select(sample, everything())

### Format TelSeq outputs
telseq$tel <- rowSums2(as.matrix(telseq[, paste0("TEL", c(6:16))]))
telseq_sum <- telseq %>%
  group_by(Library) %>%
  dplyr::summarise(Total = sum(Total),
            Mapped = sum(Mapped),
            telomere = sum(tel),
            reads_gc = sum(GC4))
telseq_sum$tel_norm <- (telseq_sum$telomere*1000000)/telseq_sum$reads_gc

### Write tables
write.table(contexts_norm, file.path(outdir, paste0(project, "_TelomereHunter_contexts_normalized.txt")), sep = "\t", row.names = FALSE)
write.table(contexts_per, file.path(outdir, paste0(project, "_TelomereHunter_contexts_percent.txt")), sep = "\t", row.names = FALSE)
write.table(summaries, file.path(outdir, paste0(project, "_TelomereHunter_summaries.txt")), sep = "\t", row.names = FALSE)
write.table(telseq_sum, file.path(outdir, paste0(project, "_TelSeq.txt")), sep = "\t", row.names = FALSE)





