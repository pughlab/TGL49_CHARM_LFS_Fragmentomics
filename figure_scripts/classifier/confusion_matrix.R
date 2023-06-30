library(ggplot2)
library(ggpubr)

outdir <- "/Users/derekwong/Library/CloudStorage/GoogleDrive-derekwong90@gmail.com/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier"

### Make confusion matrices
data <- data.frame(A = c("Negative", "Negative", "Positive", "Positive"), #Clinical
                   B = c("Negative", "Positive", "Negative", "Positive"), #ctDNA
                   Y = c(17, 2, 6, 8),
                   color = c("neg", "false_pos", "false_neg", "pos"))
data$A <- factor(data$A, levels = c("Positive", "Negative"))
data$B <- factor(data$B, levels = c("Negative", "Positive"))
totals <-data %>% 
  group_by(A) %>%
  dplyr::summarise(total = sum(Y))
data <- merge(data, totals, by = "A")
data$per <- round(data$Y/data$total, 4)*100
ppv <- format(round(data$Y[data$A == "Positive" & data$B == "Positive"]/sum(data$Y[data$B == "Positive"])*100, 2), nsmall = 2)
npv <- format(round(data$Y[data$A == "Negative" & data$B == "Negative"]/sum(data$Y[data$B == "Negative"])*100, 2), nsmall = 2)

data2 <- data.frame(A = c("Negative", "Negative", "Positive", "Positive"), #Clinical
                   B = c("Negative", "Positive", "Negative", "Positive"), #ctDNA
                   Y = c(17, 2, 4, 10),
                   color = c("neg", "false_pos", "false_neg", "pos"))
data2$A <- factor(data2$A, levels = c("Positive", "Negative"))
data2$B <- factor(data2$B, levels = c("Negative", "Positive"))
totals <-data2 %>% 
  group_by(A) %>%
  dplyr::summarise(total = sum(Y))
data2 <- merge(data2, totals, by = "A")
data2$per <- round(data2$Y/data2$total, 4)*100
ppv2 <- format(round(data2$Y[data2$A == "Positive" & data2$B == "Positive"]/sum(data2$Y[data2$B == "Positive"])*100, 2), nsmall = 2)
npv2 <- format(round(data2$Y[data2$A == "Negative" & data2$B == "Negative"]/sum(data2$Y[data2$B == "Negative"])*100, 2), nsmall = 2)

### Set theme
theme <- theme(plot.title = element_text(hjust = 0.5, size = 13), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = "none",
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               axis.text = element_text(size = 13),
               axis.title = element_text(size = 13),
               axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
               strip.background = element_blank(),
               strip.text = element_text(size = 13))

colors <- c(pos = "#B2DF8A", neg = "grey65", false_pos = "#FDBF6F", false_neg = "#FB9A99")

### Plot cancer positives
mat <- ggplot(data, aes(A, B)) +
  geom_tile(aes(fill = per), color = NA, size = 1, alpha = 0.5) +
  geom_text(aes(label = paste0(c("True\nNegative", "False\nPositive", "False\nNegative", "True\nPositive"), "\n", 
                               "(", Y, ")", "\n", per, "%")), vjust = 0.5, size = 3.5, lineheight = 1) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(0, 100)) +
  xlab(paste0("Clinical Cancer Diagnosis\n\nPPV = ", ppv, "%\n","NPV = ", npv, "%")) +
  ylab("cfDNA Fragmentation") +
  ggtitle("TP53m-carriers") +
  theme +
  coord_equal()
mat

ggsave(file.path(outdir, "confusion_matrix.pdf"), mat, width = 3.25, height = 4)

mat2 <- ggplot(data2, aes(A, B)) +
  geom_tile(aes(fill = per), color = NA, size = 1, alpha = 0.5) +
  geom_text(aes(label = paste0(c("True\nNegative", "False\nPositive", "False\nNegative", "True\nPositive"), "\n", 
                               "(", Y, ")", "\n", per, "%")), vjust = 0.5, size = 3.5, lineheight = 1) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(0, 100)) +
  xlab(paste0("Clinical Cancer Diagnosis\n\nPPV = ", ppv2, "%\n","NPV = ", npv2, "%")) +
  ylab("cfDNA Fragmentation") +
  ggtitle("TP53m-carriers") +
  theme +
  coord_equal()
mat2

ggsave(file.path(outdir, "confusion_matrix2.pdf"), mat2, width = 3.25, height = 4)
