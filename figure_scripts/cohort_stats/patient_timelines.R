library(tidyverse)
library(ggpubr)

### Set variables
path <- ""
outdir <- ""

data <- read.delim("data/patient_timelines.txt")

samples <- read.delim(file.path(path, "sample_list.txt"))
exclude <- c("TGL49_0025_Cf_U_PE_321_WG", "TGL49_0035_Cf_U_PE_310_WG", "TGL49_0041_Cf_U_PE_317_WG", "TGL49_0209_Cf_U_PE_373_WG")
samples <- samples[!(samples$sWGS %in% exclude), ]
samples <- samples[!(samples$notes == "NS"), ]
samples <- samples[!(samples$sWGS == ""), ]

data <- merge(data, samples[, c("ext_ID", "timepoint")], by = c("ext_ID", "timepoint"))

### Find latest timepoint
latest_timepoint <- data %>%
  group_by(ext_ID) %>%
  dplyr::summarise(last = max(blood),
                   n = n())
latest <- rep(latest_timepoint$last, latest_timepoint$n)
data$last <- latest

### Find patients with serial samples
data_serial <- data[data$ext_ID %in% data$ext_ID[duplicated(data$ext_ID)],]
serial_patients <- unique(data_serial$ext_ID)

data_sero <- distinct(data_serial, ext_ID, cancer_status, .keep_all= TRUE)
homo_sero <- data_sero[!(data_sero$ext_ID %in% data_sero$ext_ID[duplicated(data_sero$ext_ID)]), ]

all_negative <- homo_sero$ext_ID[homo_sero$cancer_status == "negative"]
all_positive <- homo_sero$ext_ID[homo_sero$cancer_status == "positive"]

data_sero <- data_sero[data_sero$ext_ID %in% data_sero$ext_ID[duplicated(data_sero$ext_ID)],]

data_sero <- data_sero[order(data_sero$ext_ID,
                             data_sero$cancer_status), ]

sero_multi <- c("LFS3", "LFS5", "LFS40")
data_sero <- data_sero[!(data_sero$ext_ID %in% c("LFS3", "LFS5", "LFS40")), ]

data_sero$timepoint <- as.numeric(data_sero$timepoint)

sero_neg_pos <- data_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
sero_neg_pos <- sero_neg_pos[sero_neg_pos$difference > 0, ]
forward_phenoconverters <- unique(sero_neg_pos$ext_ID)

sero_pos_neg <- data_sero %>% 
  group_by(ext_ID) %>% 
  dplyr::mutate(difference = c(diff(timepoint)))
sero_pos_neg <- sero_pos_neg[sero_pos_neg$difference < 0, ]
reverse_phenoconverters <- unique(sero_pos_neg$ext_ID)

### Set patient Types
AC <- unique(data$ext_ID[data$cancer_status == "positive"])
PC <- unique(data$ext_ID[data$cancer_status == "negative" & data$previous_cancer == "yes"])
PC <- PC[!(PC %in% AC)]
H <- unique(data$ext_ID[data$cancer_status == "negative" & data$previous_cancer == "no"])
H <- H[!(H %in% AC)]

data$type <- ifelse(data$ext_ID %in% AC, "LFS-AC",
                    ifelse(data$ext_ID %in% PC, "LFS-PC", "LFS-H"))

data$sample_type <- ifelse(data$cancer_status == "positive", "LFS-AC",
                           ifelse(data$cancer_status == "negative" & data$previous_cancer == "yes", "LFS-PC", "LFS-H"))

### Format data
data$type <- factor(data$type, levels = c("LFS-AC", "LFS-PC", "LFS-H"))
data$sample_type <- factor(data$sample_type, levels = c("LFS-AC", "LFS-PC", "LFS-H"))
data$cancer_status <- factor(data$cancer_status, levels = c("negative", "positive"),
                             labels = c("Negative", "Positive"))

data <- data[order(data$last), ]
order <- unique(data$ext_ID)
data$ext_ID <- factor(data$ext_ID, levels = order)

### Make plotting data for serial samples only
data_serial <- data[data$ext_ID %in% serial_patients, ]
data_serial$converter <- ifelse(data_serial$ext_ID %in% sero_multi, "Multiple\nCancers",
                                ifelse(data_serial$ext_ID %in% forward_phenoconverters, "Forward\nPhenoconverter",
                                       ifelse(data_serial$ext_ID %in% reverse_phenoconverters, "Reverse\nPhenoconverter",
                                              ifelse(data_serial$ext_ID %in% all_negative, "Cancer Free",
                                                     ifelse(data_serial$ext_ID %in% all_positive, "Active Cancer", NA)))))
data_serial$converter <- factor(data_serial$converter, levels = c("Multiple\nCancers", "Forward\nPhenoconverter", "Reverse\nPhenoconverter", "Active Cancer", "Cancer Free"))
data_serial <- data_serial[order(data_serial$converter,
                                 data_serial$last), ]
order <- unique(data_serial$ext_ID)
data_serial$ext_ID <- factor(data_serial$ext_ID, levels = order)

### Plot data
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               strip.background = element_rect(fill = NA),
               strip.text = element_text(size = 10),
               legend.position = c(0.65, 0.25),
               legend.key = element_blank(),
               axis.text = element_text(size = 10),
               axis.text.y = element_blank(),
               axis.title = element_text(size = 12))

plot_status <- ggplot(data, aes(blood, ext_ID, group = ext_ID, color = cancer_status)) + 
  geom_line(color = "grey") +
  geom_point(pch = 16, size = 2) +
  ggtitle("Samples by Cancer Status") +
  xlab("Months") +
  ylab("Patient") +
  labs(color = "Cancer Status") +
  scale_color_manual(values = c("black","red")) +
  theme
plot_status

plot_type <- ggplot(data, aes(blood, ext_ID, group = ext_ID, color = sample_type)) + 
  geom_line(color = "grey") +
  geom_point(pch = 16, size = 2) +
  ggtitle("Samples by Patient Status") +
  xlab("Months") +
  ylab("Patient") +
  labs(color = "Patient Status") +
  scale_color_manual(values = c("#E31A1C", "#33A02C", "#1F78B4")) +
  theme
plot_type

plot_pheno_status <- ggplot(data_serial, aes(blood, ext_ID, group = ext_ID, color = cancer_status)) + 
  geom_line(color = "grey") +
  geom_point(pch = 16, size = 2) +
  ggtitle("Cancer Status") +
  xlab("Months") +
  ylab("Patient") +
  labs(color = "Cancer Status") +
  scale_color_manual(values = c("black","red")) +
  facet_grid(converter~., scales = "free_y", space = "free_y") +
  theme +
  theme(strip.text.y = element_blank(),
        legend.position = c(0.65, 0.6))
plot_pheno_status

plot_pheno_type <- ggplot(data_serial, aes(blood, ext_ID, group = ext_ID, color = sample_type)) + 
  geom_line(color = "grey") +
  geom_point(pch = 16, size = 2) +
  ggtitle("Patient Status") +
  xlab("Months") +
  ylab("Patient") +
  labs(color = "Patient Status") +
  scale_color_manual(values = c("#E31A1C", "#33A02C", "#1F78B4")) +
  facet_grid(converter~., scales = "free_y", space = "free_y") +
  theme +
  theme(strip.text.y = element_text(angle = 0),
        legend.position = c(0.65, 0.6))
plot_pheno_type

Figure <- ggarrange(plot_status, plot_type, nrow = 1, align = "hv")
Figure

Figure2 <- ggarrange(plot_pheno_status, plot_pheno_type, nrow = 1, align = "hv")
Figure2

ggsave(file.path(outdir, "samples_timelines.pdf"), Figure, width = 7, height = 8)
ggsave(file.path(outdir, "longitudinal_timelines.pdf"), Figure2, width = 7, height = 5.5)
