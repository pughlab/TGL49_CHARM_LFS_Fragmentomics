library(tidyverse)
library(networkD3)
library(webshot)
library(htmlwidgets)
library(ComplexHeatmap)

### Set variables
path <- "/Users/derekwong/Library/CloudStorage/OneDrive-UHN/Post-Doc/CHARM_Project/LFS"
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_clinical/figures/cohort stats"

data <- read.delim(file.path(path, "samples/sample_list.txt"))
data <- data[!(data$sWGS %in% c("TGL49_0041_Cf_U_PE_317_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0209_Cf_U_PE_373_WG")), ] ### these patients excluded due to CHIP/sample swaps
data <- data[!(data$cfMeDIP %in% c("TGL49_0316_Pl_T_PE_301_CM", "TGL49_0060_Cf_n_PE_301_CM")), ] # These samples excluded (unsure origins)

### Keep serial phenoconverters samples
data <- data[!(data$notes == "NS"), ]
data <- data[data$ext_ID %in% data$ext_ID[duplicated(data$ext_ID)], ]

### positive stats
positive <- data[data$cancer_status == "positive", ]
positive <- positive[!duplicated(positive[, c("ext_ID", "cancer_type")]),]

tumors_adult <- nrow(positive[positive$Age == "adult", ])
tumors_ped <- nrow(positive[positive$Age == "pediatric", ])

pos_patients_adult <- length(unique(positive$ext_ID[positive$Age == "adult"]))
pos_patients_ped <- length(unique(positive$ext_ID[positive$Age == "pediatric"]))

pos_adult_id <- unique(positive$ext_ID[positive$Age == "adult"])
pos_ped_id <- unique(positive$ext_ID[positive$Age == "pediatric"])

### negative stats
negative <- data[data$cancer_status == "negative", ]

neg_samples_adult <- nrow(negative[negative$Age == "adult", ])
neg_samples_ped <- nrow(negative[negative$Age == "pediatric", ])

neg_patients_adult <- length(unique(negative$ext_ID[negative$Age == "adult"]))
neg_patients_ped <- length(unique(negative$ext_ID[negative$Age == "pediatric"]))

neg_samples_adult_cancer <- nrow(negative[negative$Age == "adult" & negative$ext_ID %in% pos_adult_id, ])
neg_samples_ped_cancer <- nrow(negative[negative$Age == "pediatric" & negative$ext_ID %in% pos_ped_id, ])

adults <- length(unique(negative$ext_ID[negative$Age == "adult"]))
pediatric <- length(unique(negative$ext_ID[negative$Age == "pediatric"]))
