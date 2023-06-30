library(tidyverse)
library(plyr)
library(GeneCycle)
library(matrixStats)

#analyses <- c("DHS", "TCGA", "TP53", "Immune", "uveal", "housekeeping", "TFBS", "TFBS_Ulz") ### Others = "hematopoietic"
analyses <- c("LFS")
### Set working variables
for (analysis in analyses) {
  ### If only want to parse one specific analysis
  #analysis <- "Immune_Calderon"
  path <- file.path("/Users/derekwong/Desktop/H4H/projects/CHARM/LFS/griffin/output/nucleosome_profiling", analysis, "coverage")
  outdir <- file.path("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/griffin", analysis)
  project <- "CHARM_LFS"
  
  ### Make outdir
  dir.create(outdir, showWarnings = FALSE)
  
  ### Get the site folders in the directory
  sites <- list.dirs(path, full.names = FALSE)
  discard <- c("", "stop")
  sites <- sites[!(sites %in% discard)]
  
  ### Create a loop to extract files and data from all sites
  for(site in sites){ 
    
    ### Get the list of files
    filenames <- list.files(file.path(path, site), pattern = "coverage.txt", recursive = TRUE, full.names = TRUE)
    
    ### Read in data and rbind into one dataframe
    datalist <- lapply(filenames, function(x){read.delim(file = x, header = TRUE)})
    griffin <- do.call(rbind, datalist)
    
    ### Make header for distances (they are not in order in the original files)
    header <- filenames[1]
    header <- read.delim(header, header = FALSE)
    header <- header[1, ]
    header <- header %>% dplyr::select(where(is.numeric))
    header <- t(header)
    
    ### Format griffin coverages table
    keep <- c("GC_corrected_total_read_ends", "GC_corrected_total_read_starts",
              "background_normalization", "endpoint", "number_of_sites", "sample", 
              "site_name", "total_read_ends", "total_read_starts")
    
    ### Subset the griffin analysis summaries
    summary <- griffin[griffin$GC_correction == "GC_corrected", 
                       colnames(griffin) %in% keep]
    
    ### Subset the raw uncorrected values
    raw <- griffin[griffin$GC_correction == "none", ]
    row.names(raw) <- raw$sample
    raw <- raw[ , !(colnames(raw) %in% keep)]
    raw <- raw[ , !(colnames(raw) %in% c("GC_correction"))]
    raw <- as.data.frame(t(raw))
    raw$distance <- header
    raw <- raw %>% dplyr::select(distance, everything())
    raw <- raw[order(raw$distance), ]
    
    ### Subset the GC corrected values
    corrected <- griffin[griffin$GC_correction == "GC_corrected", ]
    row.names(corrected) <- corrected$sample
    corrected <- corrected[ , !(colnames(corrected) %in% keep)]
    corrected <- corrected[ , !(colnames(corrected) %in% c("GC_correction"))]
    corrected <- as.data.frame(t(corrected))
    corrected$distance <- header
    corrected <- corrected %>% dplyr::select(distance, everything())
    corrected <- corrected[order(corrected$distance), ]
    
    ### Calculate the features for each sample (mean coverage, midpoint coverage, amplitude)
    stats <- corrected[, !(colnames(corrected) %in% "distance")]
    mean <- colMeans(stats)
    midpoint <- colMeans(stats[rownames(stats) %in% c("X.30", "X.15", "0", "X15", "X30"), ])
    fft <- GeneCycle::periodogram(stats)[["spec"]]
    fft <- colMaxs(fft)
    stats <- rbind(mean, midpoint, fft)
    features <- c(paste0(site, "_coverage"), paste0(site, "_midpoint"), paste0(site, "_amp"))
    stats <- cbind(features, stats)
    
    ### Write the tables to the output directory
    write.table(summary, file.path(outdir, paste0(project, "_griffin_summary_", site, ".txt")), row.names = FALSE, sep = "\t")
    write.table(raw, file.path(outdir, paste0(project, "_griffin_raw_", site, ".txt")), row.names = FALSE, sep = "\t")
    write.table(corrected, file.path(outdir, paste0(project, "_griffin_corrected_", site, ".txt")), row.names = FALSE, sep = "\t")
    write.table(stats, file.path(outdir, paste0(project, "_griffin_features_", site, ".txt")), row.names = TRUE, sep = "\t")
  }
}
