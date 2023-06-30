library(tidyverse)
library(plyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/hereditary/projects/CHARM/LFS/read_group/output"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/read_group"
project <- "CHARM_LFS"

### Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## find files
filenames <- list.files(path, pattern = "*.txt", full.names = TRUE)
names <- basename(filenames)
names <- sub("(_WG).*", '\\1', names)
datalist <- lapply(filenames, function(x) {
  read.delim(x, row.names = NULL, header = FALSE)})

data <- data.frame()
for (i in 1:length(datalist)) {
  x <- datalist[[i]]
  name <- names[[i]]
  
  x$sample <- name
  
  data <- rbind(data, x)
}

data$V1 <- substr(data$V1, 0, 34)

write.table(data, file.path(outdir, paste0(project, "_read_groups.txt")), row.names = FALSE, sep = "\t")


