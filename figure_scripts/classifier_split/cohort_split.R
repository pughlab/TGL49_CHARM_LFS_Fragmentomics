library(tidyverse)

### Read in sample lists and format
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

data_samples <- data_samples[!(data_samples$sWGS == ""), ]
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
data_samples <- data_samples[!(data_samples$sWGS %in% exclude), ]

healthy_samples <- read.delim("/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/Healthy_control_cohorts/CHARM_HBC/sample_list.txt")
healthy_samples <- healthy_samples[!(healthy_samples$sWGS == ""), ]

negatives <- data_samples[data_samples$cancer_status == "negative", ]
positives <- data_samples[data_samples$cancer_status == "positive", ]

### Subsample cohorts
set.seed(123)
picked <- sample(seq_len(nrow(healthy_samples)), size = floor(0.8*nrow(healthy_samples)))
hbc_test <- healthy_samples[picked,]
hbc_val <- healthy_samples[-picked,]

set.seed(123)
picked <- sample(seq_len(nrow(negatives)), size = floor(0.8*nrow(negatives)))
neg_test <- negatives[picked,]
neg_val <- negatives[-picked,]

set.seed(123)
picked <- sample(seq_len(nrow(positives)), size = floor(0.8*nrow(positives)))
pos_test <- positives[picked,]
pos_val <- positives[-picked,]

### Format sample lists
samples <- bind_rows(data.frame(sWGS = neg_test$sWGS,
                                cancer_status = "negative",
                                split = "test"),
                     data.frame(sWGS = neg_val$sWGS,
                                cancer_status = "negative",
                                split = "validation"),
                     data.frame(sWGS = pos_test$sWGS,
                                cancer_status = "positive",
                                split = "test"),
                     data.frame(sWGS = pos_val$sWGS,
                                cancer_status = "positive",
                                split = "validation"),
                     data.frame(sWGS = hbc_test$sWGS,
                                cancer_status = "healthy",
                                split = "test"),
                     data.frame(sWGS = hbc_val$sWGS,
                                cancer_status = "healthy",
                                split = "validation"))

### Save sample list
write.table(samples, file.path("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier_split", "samples_split.txt"), sep = "\t", row.names = FALSE)



