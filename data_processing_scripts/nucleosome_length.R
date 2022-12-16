library(tidyverse)
library(dplyr)
library(plyr)
library(reshape2)
library(matrixStats)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/nucleosome_peaks/output/length"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/nucleosome_peaks"
project <- "TGL49_LFS"

### Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### Read in files
names <- list.files(path = path, pattern = "*proximal.txt", full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

filenames <- list.files(path = path, pattern = "*proximal.txt", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.delim(file = x)})
data_prox <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)
colnames(data_prox) <- c("length", names)

filenames <- list.files(path = path, pattern = "*distal.txt", full.names = TRUE)
datalist <- lapply(filenames, function(x){read.delim(file = x)})
data_dist <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)
colnames(data_dist) <- c("length", names)

### Calculate fragment frequencies
length <- c(1:600)
freq_prox <- data_prox[, -1]
sums_prox <- colSums2(as.matrix(freq_prox))
freq_prox <- as.data.frame(t(t(freq_prox)/sums_prox*100))
freq_prox$length <- length
freq_prox <- freq_prox %>% dplyr::select(length, everything())
write.table(freq_prox, file.path(outdir, paste0(project, "_freq_proximal.txt")), row.names = FALSE, sep = "\t")

freq_dist <- data_dist[, -1]
sums_dist <- colSums2(as.matrix(freq_dist))
freq_dist <- as.data.frame(t(t(freq_dist)/sums_dist*100))
freq_dist$length <- length
freq_dist <- freq_dist %>% dplyr::select(length, everything())
write.table(freq_dist, file.path(outdir, paste0(project, "_freq_distal.txt")), row.names = FALSE, sep = "\t")

### Check frequency distributions
data_melt_prox <- reshape2::melt(freq_prox, id = "length")
data_melt_prox$length <- as.numeric(data_melt_prox$length)

data_melt_dist <- reshape2::melt(freq_dist, id = "length")
data_melt_dist$length <- as.numeric(data_melt_dist$length)

freq_plot1 <- ggplot(data_melt_prox, aes(x = length, y = value, group = variable)) +
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free", nrow = 6) +
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
freq_plot1

freq_plot2 <- ggplot(data_melt_dist, aes(x = length, y = value, group = variable)) +
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free", nrow = 6) +
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
freq_plot2
