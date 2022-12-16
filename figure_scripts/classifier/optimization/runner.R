### Set theme
scaleFUN <- function(x) sprintf("%.2f", x)
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
               axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.key = element_rect(fill = "white"),
               legend.text = element_text(size = 12),
               legend.position = "none",
               legend.background = element_blank(),
               axis.text = element_text(size = 12),
               axis.title = element_text(size = 12),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### Set models to run
#Features.CVparam <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = TRUE, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE)
classifiers <- c("rf", "glm", "gbm", "svmRadial", "knn")
names <- c("RF", "GLM", "GBM", "SVM", "KNN")

### Run tests (healthy vs LFS)
results <- data.frame(test = c(), kappa = c(), CI = c())

for (k in c(1:length(classifiers))){
  ### Make input dataframe
  Data <- data[row.names(data) %in% samples_neg$sWGS, ]
  Data <- bind_rows(Data, normal)
  
  y <- c(rep("positive", sum(row.names(Data) %in% samples_neg$sWGS)))
  y <- c(y, rep("negative", sum(!(row.names(Data) %in% samples_neg$sWGS))))
  y <- factor(y, levels = c("negative", "positive"))
  
  ### Set variables
  class <- classifiers[[k]]
  name <- names[[k]]
  
  ### Run suite of classifiers
  set.seed(123)
  source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/classifier_opt.R")
  
  ### Get Kappa performance
  mean <- mean(kappa)
  margin <- qt(0.975, df = length(kappa) - 1)*sd(kappa)/sqrt(length(kappa))
  
  ### Apphend to dataframe
  a <- data.frame(test = name, kappa = mean, error = margin)
  results <- rbind(results, a)
}

### Plot Kappa values
results <- results[order(results$kappa), ]
results$test <- factor(results$test, levels = results$test)
plot <- ggplot(results) +
  geom_point(aes(x = test, y = kappa)) +
  geom_errorbar(aes(x = test, ymin = kappa - error, ymax = kappa + error), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  xlab("Algorithm") + 
  ylab("Kappa") +
  labs(color = "") + 
  ggtitle(title2) + 
  theme +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip()
plot
ggsave(file.path(outdir, paste0("kappa_", title, "_status.pdf")), plot, width = 2, height = 2.5)

### Run tests (LFS cancer negative vs positive)
results <- data.frame(test = c(), kappa = c(), CI = c())

for (k in c(1:length(classifiers))){
  ### Make table
  Data <- data[row.names(data) %in% data_samples$sWGS, ]
  
  y <- data_samples$cancer_status
  y <- factor(y, levels = c("negative", "positive"))
  
  ### Set variables
  class <- classifiers[[k]]
  name <- names[[k]]
  
  ### Run suite of classifiers
  set.seed(123)
  source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/classifier_opt.R")
  
  ### Get Kappa performance
  mean <- mean(kappa)
  margin <- qt(0.975, df = length(kappa) - 1)*sd(kappa)/sqrt(length(kappa))
  
  ### Apphend to dataframe
  a <- data.frame(test = name, kappa = mean, error = margin)
  results <- rbind(results, a)
}

### Plot Kappa values
results <- results[order(results$kappa), ]
results$test <- factor(results$test, levels = results$test)
plot <- ggplot(results) +
  geom_point(aes(x = test, y = kappa)) +
  geom_errorbar(aes(x = test, ymin = kappa - error, ymax = kappa + error), width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  xlab("Algorithm") + 
  ylab("Kappa") +
  labs(color = "") + 
  ggtitle(title2) + 
  theme +
  scale_y_continuous(labels = scaleFUN) +
  coord_flip()
plot
ggsave(file.path(outdir, paste0("kappa_", title, "_lfs.pdf")), plot, width = 2, height = 2.5)

### Run tests (LFS Germline Mutation)
#results <- data.frame(test = c(), kappa = c(), CI = c())
#data_samples <- data_samples[!(data_samples$mutation_type %in% c("", "1", "5")), ]
#data_samples$mutation_type <- factor(data_samples$mutation_type, levels = c("LOF", "Splice", "2", "3"),
#                                     labels = c("LOF", "Splice", "C2", "C3"))
#
#mutations <- unique(data_samples$mutation_type)
#
#for (mutation in mutations) {
#  samples <- data_samples
#  samples$mutation_type <- ifelse(samples$mutation_type == mutation, "positive", "negative")
#  
#  Data <- data[row.names(data) %in% samples$sWGS, ]
#  
#  y <- samples$mutation_type
#  y <- factor(y, levels = c("negative", "positive"))
#  
#  for (k in c(1:length(classifiers))){
#    ### Set variables
#    class <- classifiers[[k]]
#    name <- names[[k]]
#    
#    if(mutation == "C2" & class == "gbm") next
#    ### Run suite of classifiers
#    set.seed(123)
#    source("/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/classifier/optimization/classifier_opt.R")
#    
#    ### Get Kappa performance
#    mean <- mean(kappa)
#    margin <- qt(0.975, df = length(kappa) - 1)*sd(kappa)/sqrt(length(kappa))
#    
#    ### Apphend to dataframe
#    a <- data.frame(test = name, class = mutation, kappa = mean, error = margin)
#    results <- rbind(results, a)
#  }
#}

### Plot Kappa values
#results <- results[order(results$kappa), ]
#results$test <- factor(results$test, levels = results$test)
#plot <- ggplot(results) +
#  geom_point(aes(x = test, y = kappa)) +
#  geom_errorbar(aes(x = test, ymin = kappa - error, ymax = kappa + error), width = 0) +
#  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
#  xlab("Algorithm") + 
#  ylab("Kappa") +
#  labs(color = "") + 
#  ggtitle(title2) + 
#  theme +
#  scale_y_continuous(labels = scaleFUN) +
#  coord_flip()
#plot
#ggsave(file.path(outdir, paste0("kappa_", title, "_mutation.pdf")), plot, width = 2, height = 2.5)
