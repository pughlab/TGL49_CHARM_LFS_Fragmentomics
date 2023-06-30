library(dplyr)
library(matrixStats)
library(psych)
library(ComplexHeatmap)

### Set variables
path <- ""
outdir <- ""

samples <- file.path("sample_list.txt")

### Import data
samples <- read.delim(samples)

### Remove failed samples
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]
samples <- samples[samples$mutation_type %in% c("2", "3", "LOF", "Splice"), ]
samples <- samples[!(samples$notes == "NS"), ]

samples$ped_cancer <- ifelse(samples$previous_type %like% "adrenal" |
                               samples$previous_type %like% "osteo" |
                               samples$previous_type %like% "rhabdo" |
                               samples$previous_type %like% "astro" |
                               samples$previous_type %like% "glioma" , "yes", "no")

samples$adult_cancer <- ifelse(samples$previous_type %like% "breast" |
                                 samples$previous_type %like% "sarcoma" |
                                 samples$previous_type %like% "prostate" |
                                 samples$previous_type %like% "lung" |
                                 samples$previous_type %like% "leukemia" |
                                 samples$previous_type %like% "carcinoma" |
                                 samples$previous_type %like% "thyroid" |
                                 samples$previous_type %like% "phylloides" |
                                 samples$previous_type %like% "endometrial", "yes", "no")
samples$adult_cancer <- ifelse(samples$Age == "adult", samples$adult_cancer, "no")
samples$history <- ifelse(samples$ped_cancer == "yes", "pediatric onset",
                          ifelse(samples$adult_cancer == "yes" &
                                   samples$ped_cancer == "no", "adult onset", "no cancer"))

samples_neg <- samples[samples$cancer_status == "negative", ]

patients <- samples %>%
  group_by(ext_ID) %>%
  top_n(1, timepoint)

### Summarize data based on mutation clusters (samples)
data <- samples_neg %>% 
  group_by(mutation_type) %>%
  dplyr::summarise(N = n(),
                   patients = n_distinct(ext_ID),
                   previvor = sum(previous_cancer == "no"),
                   survivor = sum(previous_cancer == "yes"),
                   pediatric = sum(Age == "pediatric"),
                   adult = sum(Age == "adult"))

chisq <- chisq.test(data[, c("survivor", "previvor")])$p.value
chisq_age <- chisq.test(data[, c("adult", "pediatric")])$p.value

data2 <- bind_rows(data[1, ], colSums(data[2:4, 2:7]))
chisq2 <- chisq.test(data2[, c("survivor", "previvor")])$p.value
chisq_age2 <- chisq.test(data2[, c("adult", "pediatric")])$p.value

### Summarize data based on mutation clusters (patients)
data_patient <- patients %>% 
  group_by(mutation_type) %>%
  dplyr::summarise(N = n(),
                   patients = n_distinct(ext_ID),
                   previvor = sum(previous_cancer == "no"),
                   survivor = sum(previous_cancer == "yes"),
                   pediatric = sum(Age == "pediatric"),
                   adult = sum(Age == "adult"))
patients$cluster <- ifelse(patients$mutation_type == "2", "Missense 2", "Cohort")
input <- table(patients$history, patients$mutation_type)
a <- chisq.test(input)$p.value

data_plot <- as.data.frame(table(patients$mutation_type, patients$history))
sums <- data_plot %>%
  group_by(Var1) %>%
  dplyr::summarise(sum = sum(Freq))
data_plot$Freq <- data_plot$Freq/sums$sum*100




