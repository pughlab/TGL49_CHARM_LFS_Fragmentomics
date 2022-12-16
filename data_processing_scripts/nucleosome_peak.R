library(tidyverse)
library(dplyr)

### Set working variables
analysis <- "CTCF"
path <- file.path("/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/nucleosome_peaks/output", analysis)
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/nucleosome_peaks"
project <- "CHARM_LFS"

filenames <- list.files(path = path, pattern = ".txt", recursive = TRUE, full.names = TRUE)
names <- list.files(path = path, pattern = ".txt", recursive = TRUE, full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in data and combine (Percentages)
datalist <- lapply(filenames, function(x) {read.delim(file = x)})
names(datalist) <- names
data <- data.frame(matrix(ncol = 2003, nrow = 0))
colnames(data) <- c("sample", "chr", -1000:1000)

for (name in names) {
  x <- datalist[[name]]
  x <- as.data.frame(t(x))
  colnames(x) <- x[1,]
  x$chr <- row.names(x)
  x$sample <- name
  x <- x %>%
    select(chr, everything())
  x <- x %>%
    select(sample, everything())
  x <- x[-1,]
  data <- rbind(data, x)
}

### Write files
write.table(data, file.path(outdir, paste0(project, "_", analysis, "_peak_distances.txt")), row.names = FALSE, sep = "\t")
