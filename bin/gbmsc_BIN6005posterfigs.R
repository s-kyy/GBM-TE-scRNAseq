#!usr/bin/env Rscript

# salloc --time=2:0:0 --ntasks=1 --mem-per-cpu=12G --account=def-ytanaka
# module load r/4.0.2
# module load gcc/9.3.0
# module load gdal/3.2.3
# srun R --no-save

## Import arguements
source("~/scratch/gete-gbm/bin/util.R")
resultsPath <- "~/scratch/gete-gbm/results/"

## Set library path
.libPaths(c("/scratch/samkyy/gete-gbm/renv/library/R-4.0/x86_64-pc-linux-gnu",
    "/tmp/RtmpJsRC8Z/renv-system-library",
    .libPaths()))
.libPaths()

## Import libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse)
# library(RCurl)
# library(scales)

# mkdirToday()
# [1] "Current working directory: /scratch/samkyy/gete-gbm"
# [1] "New working directory: /scratch/samkyy/gete-gbm/results/2022-02-10"

setwd(paste0(resultsPath,"2022-02-10"))
getwd()
########### gte UMAP seurat clusters no legend

gte <- readRDS("~/scratch/gete-gbm/results/2021-09-02/gte_gbmscIntUmap-subtypes.rds")
unique(gte$gbm_subtype)

gte$gbm_subtype <- factor(gte$gbm_subtype, levels = c('Classical',
                                                'Mesenchymal',
                                                'Proneural',
                                                'Classical-Mesenchymal',
                                                'Mesenchymal-Proneural',
                                                'Classical-Proneural',
                                                'Other'),
                                             labels = c('CL',
                                                'Mes',
                                                'PN',
                                                'CL-Mes',
                                                'Mes-PN',
                                                'CL-PN',
                                                'NA'))
levels(gte$gbm_subtype)

# Colours: same order as factored gbm_subtype levels. 
retroCols <- c("indianred2", "darkolivegreen3", "steelblue2", "pink1", "cyan3") # missing MES-PN and NA


# Set Identities
Idents(gte) <- "integrated_snn_res.0.3"

# Produce Plots
p1 <- DimPlot(gte, reduction = "umap", 
                group.by = "gbm_subtype", 
                cols = retroCols, 
                label = FALSE, 
                repel = TRUE) + 
            ggtitle("") + NoLegend()
p1

ggsave("gbmsc_PC20_r03_gte_subtypes_nolabels.tiff", units="in", width=8, height=8, dpi=300, compression = 'lzw')

rm(gte)
rm(retroCols)
rm(p1)

################### ge no legends

ge <- readRDS("~/scratch/gete-gbm/results/2021-09-02/ge_gbmscIntUmap-subtypes.rds")
unique(ge$gbm_subtype)

ge$gbm_subtype <- factor(ge$gbm_subtype, levels = c('Classical',
                                                'Mesenchymal',
                                                'Proneural',
                                                'Classical-Mesenchymal',
                                                'Mesenchymal-Proneural',
                                                'Classical-Proneural',
                                                'Other'),
                                             labels = c('CL',
                                                'Mes',
                                                'PN',
                                                'CL-Mes',
                                                'Mes-PN',
                                                'CL-PN',
                                                'NA'))
levels(ge$gbm_subtype)

# Colours: same order as factored gbm_subtype levels. 
colours <- c("indianred2", "darkolivegreen3", "steelblue2", "pink1", "aquamarine2", "cyan3", "gray70")


# Set Identities
Idents(ge) <- "integrated_snn_res.0.3"

# Produce Plots
p1 <- DimPlot(ge, reduction = "umap", 
                group.by = "gbm_subtype", 
                cols = colours, 
                label = FALSE, 
                repel = TRUE) + 
            ggtitle("") + NoLegend()
p1

ggsave("gbmsc_PC20_r03_ge_subtypes_nolabels.tiff", units="in", width=8, height=8, dpi=300, compression = 'lzw')

rm(ge)
rm(colours)
rm(p1)

gc(verbose = TRUE,reset=TRUE)

########## teRatio UMAP no legend

gte <- readRDS("~/scratch/gete-gbm/results/2021-09-02/gte_gbmscIntUmap-subtypes.rds")

FeaturePlot(gte, features = "teRatio", reduction = "umap", min.cutoff=19.26, label = FALSE) + ggtitle("") 
ggsave("gbmsc_PC20_r03_gte_teRatio_nolabel.tiff", units="in", width=10, height=10, dpi=300, compression = 'lzw')

gc(verbose = TRUE,reset=TRUE)

########## teRatio Violin plot with colours

df <- gte@meta.data %>% dplyr::select(c(integrated_snn_res.0.3, gbm_subtype, teRatio))

#factor gbm_subtype
df$gbm_subtype <- factor(gte$gbm_subtype, levels = c('Classical',
                                                'Mesenchymal',
                                                'Proneural',
                                                'Classical-Mesenchymal',
                                                'Classical-Proneural'),
                                             labels = c('CL', 'Mes','PN','CL-Mes', 'CL-PN'))

retroCols <- c("indianred2", "darkolivegreen3", "steelblue2", "pink1", "cyan3") # missing MES-PN and NA

# Basic violin plot (ensuring adjusted p-values are available)
p <- ggplot(df, aes(x=integrated_snn_res.0.3, y=teRatio, fill=gbm_subtype)) +
        geom_violin(trim=FALSE) + geom_boxplot(width=0.1, show.legend=FALSE, fill="white") +
        scale_color_manual(values = retroCols) +
        xlab("Cluster ID") + ylab("% Transposon Counts") + labs(fill = "GBM Subtype") + theme_classic() +
        geom_hline(yintercept = mean(df$teRatio), linetype = 2)

size <- 5

pdf("gbmsc_violin_gte_teRatio_nopval.pdf", width = size*2, height = size)
p
dev.off()

gc()

p
ggsave("gbmsc_violin_gte_teRatio_nopval.tiff", units="in", width=10, height=5, dpi=300, compression = 'lzw')

gc()
rm(p)

########## Volcano Plot

c7 <- readRDS("/home/samkyy/scratch/gete-gbm/results/2021-10-14/cluster7-TEmarkers.rds")

## Set differetially expressed genes
## add a column of NAs
c7$diffexp <- "NO"
## Set as "UP" 
c7$diffexp[c7$avg_log2FC >= 0.5 & c7$p_adj_FDR < 0.05] <- "UP"
## Set as "DOWN"
c7$diffexp[c7$avg_log2FC <= -0.5 & c7$p_adj_FDR < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
c7$delabel <- NA
c7$delabel[c7$diffexp != "NO"] <- c7$gene_id[c7$diffexp != "NO"]

length(c7$diffexp[c7$avg_log2FC >= 0.5 & c7$p_adj_FDR < 0.05]) #112

gc()

head(c7)

library(ggrepel)
options(ggrepel.max.overlaps = 15)

p <- c7 %>% ggplot(aes(x=avg_log2FC, y=nlog10p.FDR, col=diffexp, label=delabel)) +
        geom_point() + xlim(-1.2, 1.2) +
        theme_classic() + 
        geom_text_repel() + xlab("Average log2FC") + ylab("-Log10(FDR)") + labs(col = "") +
        scale_color_manual(values=c("black", "red")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="black") +
        geom_hline(yintercept=-log10(0.05), col="black") + NoLegend()

# Save figures
# size <- 5
# pdf("r_GBMSC_DEG-TE-C7.pdf", width = size, height = size)
# p
# dev.off()

ggsave("r_GBMSC_DEG-TE-C7.tiff", plot = p, units="in", width=7, height=7, dpi=300, compression = 'lzw')

rm(p)
rm(c7)

gc()

########## Gene Ontology

source("~/scratch/gete-gbm/bin/util.R")
source("~/scratch/gete-gbm/bin/util_go.R")
resultsPath <- "~/scratch/gete-gbm/results/"

## Set library path
.libPaths(c("/scratch/samkyy/gete-gbm/renv/library/R-4.0/x86_64-pc-linux-gnu",
    "/tmp/RtmpJsRC8Z/renv-system-library",
    .libPaths()))
.libPaths()

## Import libraries
library(tidyverse)

setwd(paste0(resultsPath,"2022-02-10"))
getwd()

## process cluster 7 markers from 2021-06-21
load(paste0(resultsPath,"2021-06-21/gte_markers.RData"))
head(GBM.GO)

length(gte_markers$gene[gte_markers$cluster == '7']) #8

head(gte_markers$gene[gte_markers$cluster == '7' & gte_markers$avg_log2FC >= 0.5 & gte_markers$p_val_adj < 0.05])
length(gte_markers$gene[gte_markers$cluster == '7' & gte_markers$avg_log2FC >= 0.5 & gte_markers$p_val_adj < 0.05]) #8

head(gte_markers[gte_markers$cluster == '7',])

# Saved to 2022-02-10
write.csv(gte_markers$gene[gte_markers$cluster == '7' & gte_markers$avg_log2FC >= 0.5 & gte_markers$p_val_adj < 0.05],
            "gte_markersc72021-06-21.csv", row.names= FALSE) 

rm(gte_markers)
gc()

## process cluster 7 markers DEGs from 2021-10-14
c7deg <- readRDS(paste0(resultsPath,"2021-10-14/cluster7.markers.rds"))
View(c7deg)

View(as.data.frame(rownames(c7deg[c7deg$avg_log2FC >= 0.5 & c7deg$p_adj_FDR < 0.05,])))
write.csv(as.data.frame(rownames(c7deg[c7deg$avg_log2FC >= 0.5 & c7deg$p_adj_FDR < 0.05,])),
            "gte_markersc72021-10-14_up.csv", row.names= FALSE)

## process GO terms from 2021-10-14 DEG c7 markers
GBM.GO <- read.csv(paste0(resultsPath,"2022-01-24/gbmsc_c7_GOsummary_GBM.csv"))

## all reference genes
ann <- "org.Hs.eg.db"

ref2eg <- org.Hs.egSYMBOL2EG # Gene Symbols to Entrez IDs
mapped_seqs <- mappedkeys(ref2eg) # Gene Symbols

universe <- ls(org.Hs.egSYMBOL) # list of all objects names for the human SYMBOL submap

# convert to gene clusters
gene_id <- unique(as.character(rownames(c7deg[c7deg$avg_log2FC >= 0.5 & c7deg$p_adj_FDR < 0.05,])))
entrez <- intersect(gene_id, mapped_seqs)








#### Extract clusters
cluster_id <- unique(as.character(cluster.markers$cluster))
print(cluster_id)
if(length(cluster_id) == 0){
    print("no clusters")
    return()
    gene_id <- unique(as.character(cluster.markers$gene))
    if(length(gene_id) == 0){
        print("no genes")
        return()
}
} else { 
    #### Extract genes by cluster and convert to entrez_id's
    print("extracting gene symbols by cluster")
    
    
    id <- vector("list", length = length(cluster_id))
    entrez_id <- vector("list", length = length(cluster_id))
    nonunique_id <- vector("list", length = length(cluster_id))
    
    print(paste0("Length of cluster_id : ", length(cluster_id)))

    for(i in 1:length(cluster_id)){
        n = as.numeric(cluster_id[i])
        id[[i]] <- cluster.markers$gene[which(cluster.markers$cluster == n)]
        id[[i]] <- intersect(id[[i]], mapped_seqs)

        entrez_id[[i]] <- na.omit(unique(as.list(unlist( ref2eg[id[[i]]] ) )))
        nonunique_id[[i]] <- length(id[[i]]) - length(entrez_id[[i]])
        print(paste0("Number of non-unique genes in cluster", n, " = ", nonunique_id[[i]]))
    }
    print("converted gene symbols to entrez ids")
}