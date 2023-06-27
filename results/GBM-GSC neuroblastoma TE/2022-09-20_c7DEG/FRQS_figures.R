
#deg <- read.csv("S:\\My Drive\\20##\\GRAD SCHOOL\\Grad Project Documents\\TE Experiment Data\\GBMSC\\filt_GE_clustermarkers.csv")
#View(deg)


# 2022-09-20: Figures for FRQS paper 
library(ggplot2)
library(dplyr)
library(ggrepel)
options(ggrepel.max.overlaps = 100)


set.seed(10)
setwd(choose.dir()) 
getwd() # "D:/backup files/getegbm2022-02-12/gete-gbm/results/GBM-GSC neuroblastoma TE/2022-09-20"

# load dataset
# gte7 <- readRDS("D:\\backup files\\2021-12-20\\Cluster\\scratch\\gete-gbm\\results\\GBM-GSC neuroblastomaTE\\2021-10-14\\cluster7.markers.rds") 
# View(gte7)


# Add column to gte7 to colour the dots based on type of gene (upregulated). 
# gte7 <- gte7 %>% tibble::rownames_to_column(var = "gene_id")
# head(gte7)

# te_genes <- read.table("D:\\backup files\\2021-12-20\\Cluster\\scratch\\gete-gbm\\data\\geneName_TE_v23-1gtf.txt",
#               sep=",", header = TRUE,
#               col.names=c("gene_id", "TE_class", "TE_family"), 
#               fill=FALSE, 
#               strip.white=TRUE)
# 
# te_genes$gene_id <- sub("_", "-", te_genes$gene_id) # in seurat object, `_` were replaced with `-`.
# matched_TEs <- intersect(te_genes$gene_id, gte7$gene_id)
# length(matched_TEs)
# 
# gte7 <- merge(x = gte7, y = te_genes, by = "gene_id", all.x = TRUE)
#   
# # Add col of NAs
#   gte7$deg <- "NO"
# # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
#   gte7$deg[gte7$avg_log2FC > 0.5 & gte7$p_val_adj < 0.05] <- "UP"
# # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
#   gte7$deg[gte7$avg_log2FC < -0.5 & gte7$p_val_adj < 0.05] <- "DOWN"
# 
# # Set -nlogp-values to 320 as max
#   gte7$nlog10p.adj_max <- gte7$nlog10p.adj
#   gte7$nlog10p.adj_max[gte7$nlog10p.adj == Inf] <- 320
# 
# # Plot
# windows(); ggplot(gte7, aes(x=avg_log2FC, y = nlog10p.adj_max)) + 
#   geom_point(aes(color = deg), size=2.5, show.legend = FALSE) + 
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey60")) + 
#   xlim(-3.5, 2) + ylim(-2, 325) +
#   geom_hline(yintercept = -log10(0.05), col = "darkgrey", linetype = 'dashed', size = 1.05) + 
#   geom_vline(xintercept = c(-0.5, 0.5), col = "darkgrey", linetype = 'dashed', size = 1.05) + 
#   scale_color_manual(values=c("dodgerblue", "grey", "palevioletred1")) 

# Add Labels column (manual)
  # gte7$deglabs <- NA 
  # write.csv(gte7,"gte_c7_markers_volcano.csv", row.names= FALSE)
  # 
  # gte7 %>% 
  #   filter(gte7$TE_family == "Alu" & gte7$deg == "UP") %>% 
  #   arrange(desc(nlog10p.adj_max)) %>% 
  #   select(gene_id, nlog10p.adj_max, TE_class)
  # 
  # gte7 %>% 
  #   filter(gte7$TE_family == "L1" & gte7$deg == "UP") %>% 
  #   arrange(desc(nlog10p.adj_max)) %>% 
  #   select(gene_id, nlog10p.adj_max, TE_class)
  # 
  # gte7 %>% 
  #   filter(gte7$TE_family == "SVA" & gte7$deg == "UP") %>% 
  #   arrange(desc(nlog10p.adj_max)) %>% 
  #   select(gene_id, nlog10p.adj_max, TE_class)
  
gte7 <- read.csv("gte_c7_markers_volcano.csv")

# Plot
# windows(); ggplot(gte7, aes(x=avg_log2FC, y = nlog10p.adj_max, label= deglabs)) + 

p1 <- ggplot(gte7, aes(x=avg_log2FC, y = nlog10p.adj_max, label= deglabs)) +
        geom_point(aes(color = deg), size=2.4, show.legend = FALSE) + 
        theme(panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "grey50")) + 
        xlim(-3.5, 2.5) + ylim(-2, 330) +
        geom_hline(yintercept = -log10(0.05), col = "darkgrey", linetype = 'dashed', size = 1.05) + 
        geom_vline(xintercept = c(-0.5, 0.5), col = "darkgrey", linetype = 'dashed', size = 1.05) + 
        scale_color_manual(values=c("dodgerblue", "grey", "palevioletred1")) + 
        geom_text_repel(box.padding = 1, size=4.5, xlim = c(0.50, 2.5)) + 
        xlab("") + ylab("")

p1
ggsave("gte_c7_volcano.tiff", units = "in", width=6, height=6, dpi=300, compression = 'lzw')

####################

cancer.GO <- read.csv("gte_GO_GBM.NSC.c7_filt.csv")
cancer.GO$GeneRatio <- round(cancer.GO$Count/cancer.GO$Size*100, 3)
cancer.GO <- arrange(cancer.GO, group, GeneRatio) 
cancer.GO$group <- factor(cancer.GO$group, c("TE-enriched", "GBM", "NSC"))
cancer.GO$Term <- factor(cancer.GO$Term, unique(cancer.GO$Term))
View(cancer.GO)
glimpse(cancer.GO)

p1 <- cancer.GO %>% ggplot() + geom_point(aes(x=group, y=Term, color = `nlog10p.adj`, size = GeneRatio)) + 
  theme(panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = "grey50")) +
  xlab("") + ylab("") +
  scale_color_gradient2(midpoint=8, low="skyblue", mid = "darkblue", high="darkblue")
p1 + theme(axis.text.x=element_text(angle = 45, hjust = 1))

