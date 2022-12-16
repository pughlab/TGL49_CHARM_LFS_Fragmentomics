library(tidyverse)
library(matrixStats)
library(reshape2)
library(plyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/peak_reads/output/"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/peak_reads"
project <- "CHARM_LFS"

### Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### List sites
sites <- list.dirs(path, full.names = FALSE)
sites <- sites[!(sites == "")]

for (site in sites) {
  ### List files
  filenames <- list.files(path = file.path(path, site), pattern = "*picard.txt", full.names = TRUE)
  names <- list.files(path = file.path(path, site), pattern = "*.txt", full.names = FALSE)
  names <- sub("(_WG).*", '\\1', names)
  
  ### Read in and merge
  datalist <- lapply(filenames, function(x){read.delim(file = x, skip = 13, colClasses = c(rep("integer", 2), rep("NULL", 2)), 
                                                       header = TRUE, col.names = c("length", "freq", "X", "Y"))})
  data <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)
  
  data[is.na(data)] <- 0
  
  ### Format fragment data and save
  row.names(data) <- data$length
  data <- data[,-1]
  colnames(data) <- names
  lengths <- c(10:600)
  data <- data[row.names(data) %in% lengths,]
  data$length <- row.names(data)
  data <- data %>% dplyr::select(length, everything())
  write.table(data, file.path(outdir, paste0(project, "_peak_fragments_", site ,".txt")), row.names = FALSE, sep = "\t")
  
  ### Calculate fragment frequencies
  freq <- data[, -1]
  sums <- colSums2(as.matrix(freq[row.names(freq) %in% c(10:250), ]))
  freq <- as.data.frame(t(t(freq)/sums*100))
  freq$length <- row.names(freq)
  freq <- freq %>% dplyr::select(length, everything())
  write.table(freq, file.path(outdir, paste0(project, "_peak_fragment_freq_", site, ".txt")), row.names = FALSE, sep = "\t")
  
  ### Calculate fragment proportions
  short <- freq[ , -1]
  short <- short[rownames(short) %in% c(20:150), ]
  prop <- colSums2(as.matrix(short))/100
  prop <- as.data.frame(prop)
  
  ### Format fragment proportions and save
  prop$sample <- colnames(short)
  write.table(prop, file.path(outdir, paste0(project, "_peak_proportions_", site, ".txt")), row.names = FALSE, sep = "\t")
  
  ### Check frequency distributions
  data_melt <- melt(freq, value = "length")
  data_melt$length <- as.numeric(data_melt$length)
  
  freq_plot <- ggplot(data_melt, aes(x = length, y = value, group = variable)) +
    geom_line(size = 1) +
    facet_wrap(~variable, scales = "free", ncol = 5) +
    xlab("Fragment Size") + 
    ylab("Frequency (%)") +
    ggtitle("") + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15))
  
  ggsave(file.path(outdir, paste0("Frequency_check_", site, ".pdf")), freq_plot, device = "pdf", width = 12, height = 100, units = "in", limitsize = FALSE)
}
