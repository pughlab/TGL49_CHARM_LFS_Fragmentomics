library(pacman)
library(DOSE)
library(GO.db)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(rWikiPathways)
library(RCy3)

### Download and launch Cytoscape (https://cytoscape.org/download.html)
cytoscapeVersionInfo()
cytoscapePing()

### Set paths
outdir <- "/Users/derekwong/My Drive/Post-Doc/CHARM/LFS/LFS_fragment/figures/fragment_nucleosome_TFBS"

data <- read.delim(file.path(outdir, "fragment_nucleosome_TFBS_differential.txt"))
data_lfs <- read.delim(file.path(outdir, "fragment_nucleosome_TFBS_delist.txt"))
data_samples <- read.delim("/Users/derekwong/OneDrive - UHN/Post-Doc/CHARM_Project/LFS/samples/sample_list.txt")

### Seperate out features
data$site <- gsub("_.*","", data$feature)
data$variable <- gsub(".*_","", data$feature)

### Format tables
data <- data[data$variable == "midpoint" &
               data$padj_cohort <= 0.01 &
               data$cohort < -1, 
             c("site", "cohort", "padj_cohort")]

### Convert to ENTREZ IDs
dn.genes <- data$site[data$cohort < 0]
dn.genes.entrez <- clusterProfiler::bitr(dn.genes, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

egobp <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)

head(egobp,10)
barplot(egobp, showCategory = 20)
dotplot(egobp, showCategory = 20)

egobp <- as.data.frame(egobp)
egobp_short <- egobp[1:20, ]

### Set y-axis function
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

### Set Theme
theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
               axis.title = element_text(size = 12),
               axis.line = element_line(colour = "black"),
               axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
               axis.text = element_text(size = 10),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.position = c(0.6, 0.15),
               legend.direction = "horizontal",
               legend.background = element_rect(fill = "white"),
               legend.text = element_text(size = 10),
               legend.key = element_rect(fill = "white"),
               legend.spacing.y = unit(0, "mm"),
               legend.key.size = unit(3, "mm"),
               strip.background = element_rect(color = NA, fill = NA),
               strip.text = element_text(face = "bold", size = 10))

### Plot and save
plot <- ggplot(egobp_short, aes(x = p.adjust, y = reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count, fill = Count), color = "white", pch = 21) +
  scale_fill_gradient(low = "blue", high = "red",
                      limits = c(0, 15), breaks = c(0, 5, 10, 15),
                      guide = guide_colorbar(title.position = "top",
                                             barwidth = 4,
                                             barheight = 1,
                                             label.vjust = -1,
                                             title.vjust = 1,
                                             title.hjust = 0.5)) +
  scale_size_continuous(range = c(2, 5)) +
  labs(x = "p-adjusted", y = "", fill = "# of TFs") +
  ggtitle("LFS vs Healthy") +
  guides(size = "none",) +  
  theme + 
  scale_x_continuous(trans=reverselog_trans(10))
plot
ggsave(file.path(outdir, paste0("fragment_nucleosome_TFBS_pathway_cohort.pdf")), plot, device = "pdf", width = 5.5, height = 4, units = "in")

### Pathway analysis for within LFS patients
### Convert to ENTREZ IDs
up.genes <- data_lfs$site[data_lfs$direction == "Increase"]
up.genes.entrez <- clusterProfiler::bitr(up.genes, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

dn.genes <- data_lfs$site[data_lfs$direction == "Decrease"]
dn.genes.entrez <- clusterProfiler::bitr(dn.genes, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

egobp_up <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)
head(egobp_up,10)

egobp_dn <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)
head(egobp_dn,10)

### Format plotting dataframe
egobp_up <- as.data.frame(egobp_up)
egobp_up_short <- egobp_up[1:10, ]
egobp_up_short$direction <- "Decreased Accessibility"

egobp_dn <- as.data.frame(egobp_dn)
egobp_dn_short <- egobp_dn[1:10, ]
egobp_dn_short$direction <- "Increased Accessibility"

### Plot and save
plot_lfs_dn <- ggplot(egobp_dn_short, aes(x = p.adjust, y = reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count, fill = Count), color = "white", pch = 21) +
  scale_fill_gradient(low = "blue", high = "red",
                      limits = c(0, 15), breaks = c(0, 5, 10, 15),
                      guide = guide_colorbar(title.position = "top",
                                             barwidth = 4,
                                             barheight = 1,
                                             label.vjust = -1,
                                             title.vjust = 1,
                                             title.hjust = 0.5)) +
  scale_size_continuous(range = c(2, 5)) +
  labs(x = "p-adjusted", y = "", fill = "# of TFs") +
  ggtitle("Increased Accessibility") +
  guides(size = "none",) +  
  theme + 
  theme(legend.position = c(0.75, 0.3)) +
  scale_x_continuous(trans=reverselog_trans(10))
plot_lfs_dn

plot_lfs_up <- ggplot(egobp_up_short, aes(x = p.adjust, y = reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count, fill = Count), color = "white", pch = 21) +
  scale_fill_gradient(low = "blue", high = "red",
                      limits = c(0, 15), breaks = c(0, 5, 10, 15),
                      guide = guide_colorbar(title.position = "top",
                                             barwidth = 4,
                                             barheight = 1,
                                             label.vjust = -1,
                                             title.vjust = 1,
                                             title.hjust = 0.5)) +
  scale_size_continuous(range = c(2, 5)) +
  labs(x = "p-adjusted", y = "", fill = "# of TFs") +
  ggtitle("Decreased Accessibility") +
  guides(size = "none",) +  
  theme + 
  theme(legend.position = c(0.75, 0.3)) +
  scale_x_continuous(trans=reverselog_trans(10))
plot_lfs_up

plot_lfs <- ggarrange(plot_lfs_dn, plot_lfs_up, nrow = 2, align = "hv")
plot_lfs <- annotate_figure(plot_lfs, top = text_grob("LFS Cancer Positive vs Negative", size = 16))
plot_lfs
ggsave(file.path(outdir, paste0("fragment_nucleosome_TFBS_pathway_lfs.pdf")), plot_lfs, device = "pdf", width = 7, height = 7, units = "in")
