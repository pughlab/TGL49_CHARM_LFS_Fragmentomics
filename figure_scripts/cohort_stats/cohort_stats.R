library(tidyverse)
library(networkD3)
library(webshot)
library(htmlwidgets)
library(ComplexHeatmap)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"

data_samples <- read.delim(file.path(path, "samples/sample_list.txt"))
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]
data_samples <- data_samples[!(data_samples$notes == "NS"), ]

### Find Number of patients
samples <- nrow(data_samples)

patients <- length(unique(data_samples$ext_ID))

data_pos <- data_samples[data_samples$cancer_status == "positive", ]
positive_samples <- nrow(data_pos)

data_neg <- data_samples[data_samples$cancer_status == "negative", ]
negative_samples <- nrow(data_neg)

data_previvor <- data_neg[data_neg$previous_cancer == "no", ]
data_previvor <- nrow(data_previvor)

data_survivor <- data_neg[data_neg$previous_cancer == "yes", ]
data_survivor <- nrow(data_survivor)

mut_lof <- nrow(data_samples[data_samples$mutation_type == "LOF", ])
mut_sp <- nrow(data_samples[data_samples$mutation_type == "Splice", ])
mut_1 <- nrow(data_samples[data_samples$mutation_type == "1", ])
mut_2 <- nrow(data_samples[data_samples$mutation_type == "2", ])
mut_3 <- nrow(data_samples[data_samples$mutation_type == "3", ])
mut_4 <- nrow(data_samples[data_samples$mutation_type == "4", ])
mut_5 <- nrow(data_samples[data_samples$mutation_type == "5", ])

