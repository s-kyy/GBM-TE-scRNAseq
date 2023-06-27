#module load R/3.5.0-foss-2016b-avx2
#Author: Dr. Yoshiaki Tanaka
library(Seurat)
library(Matrix)

#10x combine
combined <- Read10X("merge_brain/outs/filtered_feature_bc_matrix")
#combined_v2 <- Read10X("merge_brain_rmsk/outs/filtered_feature_bc_matrix")
#cellnames <- intersect(colnames(combined),colnames(combined_v2))
#combined <- rbind(combined[,cellnames],combined_v2[,cellnames])
#rm(combined_v2)

#seurat
combined_brain <- CreateSeuratObject(counts=combined,project="combined",names.field = 2,names.delim = "-")
orig_id <- combined_brain@meta.data$orig.ident
names(orig_id) <- rownames(combined_brain@meta.data)


combined_brain.list <- SplitObject(object=combined_brain,split.by="orig.ident")


for(i in 1:length(combined_brain.list)){
combined_brain.list[[i]] <- NormalizeData(combined_brain.list[[i]],verbose=F)
combined_brain.list[[i]] <- FindVariableFeatures(combined_brain.list[[i]],selection.method="vst",nfeatures=2500,verbose=F)
}


combined_brain.anchors <- FindIntegrationAnchors(combined_brain.list,dims=1:20)
combined_brain.intergrated <- IntegrateData(anchorset = combined_brain.anchors,dims=1:20)

library(ggplot2)
library(cowplot)

DefaultAssay(object = combined_brain.intergrated) <- "integrated"
combined_brain.intergrated <- ScaleData(object = combined_brain.intergrated, verbose = FALSE)


#library(ggplot2)
#library(cowplot)

#DefaultAssay(object = combined_brain.intergrated) <- "integrated"

combined_brain.intergrated <- RunPCA(object = combined_brain.intergrated, npcs = 20, verbose = FALSE)
combined_brain.intergrated <- RunUMAP(object = combined_brain.intergrated, reduction = "pca", dims = 1:20)
save(combined_brain.intergrated,file="combined_integrated_brain.dat")

combined_brain.intergrated <- FindNeighbors(combined_brain.intergrated,dims=1:20,reduction="pca")
combined_brain.intergrated <- FindClusters(combined_brain.intergrated)
save(combined_brain.intergrated,file="combined_integrated_brain.dat")

library(genefilter)

preclust <- combined_brain.intergrated@active.ident
clust_num <- length(unique(preclust))
ratio <- vector("list",clust_num)
names(ratio) <- 1:clust_num
ratio -> pval
ratio -> dif_gene_1p25_pval005_brain
row_count <- nrow(combined_brain.intergrated@assays$RNA)
col_count <- ncol(combined_brain.intergrated@assays$RNA)

exp <- as.matrix(combined_brain.intergrated@assays$RNA[1:5000,1:col_count])
for(i in (1:clust_num)-1){
      samp <- numeric(col_count)
      samp[which(preclust==i)] <- 1
      ratio[[(i+1)]] <- c(ratio[[(i+1)]],rowMeans(exp[,which(samp==1)]) - rowMeans(exp[,which(samp==0)]))
      pval[[(i+1)]]  <- rbind(pval[[(i+1)]],rowttests(as.matrix(exp),fac=factor(samp)))
}

exp <- as.matrix(combined_brain.intergrated@assays$RNA[5001:10000,1:col_count])
for(i in (1:clust_num)-1){
      samp <- numeric(col_count)
      samp[which(preclust==i)] <- 1
      ratio[[(i+1)]] <- c(ratio[[(i+1)]],rowMeans(exp[,which(samp==1)]) - rowMeans(exp[,which(samp==0)]))
      pval[[(i+1)]]  <- rbind(pval[[(i+1)]],rowttests(as.matrix(exp),fac=factor(samp)))
}

exp <- as.matrix(combined_brain.intergrated@assays$RNA[10001:15000,1:col_count])
for(i in (1:clust_num)-1){
      samp <- numeric(col_count)
      samp[which(preclust==i)] <- 1
      ratio[[(i+1)]] <- c(ratio[[(i+1)]],rowMeans(exp[,which(samp==1)]) - rowMeans(exp[,which(samp==0)]))
      pval[[(i+1)]]  <- rbind(pval[[(i+1)]],rowttests(as.matrix(exp),fac=factor(samp)))
}

exp <- as.matrix(combined_brain.intergrated@assays$RNA[15001:20000,1:col_count])
for(i in (1:clust_num)-1){
      samp <- numeric(col_count)
      samp[which(preclust==i)] <- 1
      ratio[[(i+1)]] <- c(ratio[[(i+1)]],rowMeans(exp[,which(samp==1)]) - rowMeans(exp[,which(samp==0)]))
      pval[[(i+1)]]  <- rbind(pval[[(i+1)]],rowttests(as.matrix(exp),fac=factor(samp)))
}

exp <- as.matrix(combined_brain.intergrated@assays$RNA[20001:25000,1:col_count])
for(i in (1:clust_num)-1){
      samp <- numeric(col_count)
      samp[which(preclust==i)] <- 1
      ratio[[(i+1)]] <- c(ratio[[(i+1)]],rowMeans(exp[,which(samp==1)]) - rowMeans(exp[,which(samp==0)]))
      pval[[(i+1)]]  <- rbind(pval[[(i+1)]],rowttests(as.matrix(exp),fac=factor(samp)))
}

exp <- as.matrix(combined_brain.intergrated@assays$RNA[25001:30000,1:col_count])
for(i in (1:clust_num)-1){
      samp <- numeric(col_count)
      samp[which(preclust==i)] <- 1
      ratio[[(i+1)]] <- c(ratio[[(i+1)]],rowMeans(exp[,which(samp==1)]) - rowMeans(exp[,which(samp==0)]))
      pval[[(i+1)]]  <- rbind(pval[[(i+1)]],rowttests(as.matrix(exp),fac=factor(samp)))
}

exp <- as.matrix(combined_brain.intergrated@assays$RNA[30001:row_count,1:col_count])
for(i in (1:clust_num)-1){
      samp <- numeric(col_count)
      samp[which(preclust==i)] <- 1
      ratio[[(i+1)]] <- c(ratio[[(i+1)]],rowMeans(exp[,which(samp==1)]) - rowMeans(exp[,which(samp==0)]))
      pval[[(i+1)]]  <- rbind(pval[[(i+1)]],rowttests(as.matrix(exp),fac=factor(samp)))
}

save(ratio,file="ratio_brain.dat")
save(pval,file="pval_brain.dat")

for(i in (1:clust_num)-1){
      dif_gene_1p25_pval005_brain[[(i+1)]] <- rownames(exp)[which(ratio[[(i+1)]] > log2(1.25) & pval[[(i+1)]][,3] < 0.05)]
}

save(dif_gene_1p25_pval005_brain,file="dif_gene_1p25_pval005_brain.dat")

rm(exp)
