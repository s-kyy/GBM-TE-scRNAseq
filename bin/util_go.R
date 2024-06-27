#!usr/bin/env Rscript

### Load Libraries

# install.packages('BiocManager')
# BiocManager::install("GOstats")
# BiocManager::install("KEGGREST")
# BiocManager::install("org.Hs.eg.db")

library(GOstats)
library(KEGGREST)
library(org.Hs.eg.db)
library(GO.db)
library(dplyr)
conflict_prefer("select", "dplyr") ## required in %>% dplyr

### Store Gene Ontology Classes and Entrez ID conversion Database
go.classes <- c("BP","MF","CC","KEGG")
hs.annotations <- "org.Hs.eg.db"
no.classes <- length(go.classes)

# Gene Symbols to Entrez IDs
ref2eg <- org.Hs.egSYMBOL2EG 
# Gene Symbols
mapped_seqs <- mappedkeys(ref2eg) 
# list of objects names for the human SYMBOL submap
universe <- ls(org.Hs.egSYMBOL) 

### Perform Gene Ontology Analysis
### cluster.markers -- a data.frame that must include $cluster and $gene columns to indicate the cluster ID and the rowname (gene symbol)

geneontology <- function(cluster.markers){
    
    #### Extract clusters
    cluster_id <- unique(as.character(cluster.markers$cluster))
    writeLines(cluster_id)
    if(length(cluster_id) == 0){
        writeLines("no clusters")
        return()
        gene_id <- unique(as.character(cluster.markers$gene))
        if(length(gene_id) == 0){
            writeLines("no genes")
            return()
    }
    } else { 
        #### Extract genes by cluster and convert to entrez_id's
        writeLines("extracting gene symbols by cluster")
        writeLines(paste0("Length of cluster_id : ", length(cluster_id)))
        
        id <- vector("list", length = length(cluster_id))
        entrez_id <- vector("list", length = length(cluster_id))
        nonunique_id <- vector("list", length = length(cluster_id))
        
        writeLines("Adding genes to each list by cluster_id")
        for(i in 1:length(cluster_id)){
            n = as.numeric(cluster_id[i])
            id[[i]] <- cluster.markers$gene[which(cluster.markers$cluster == n)]
            id[[i]] <- intersect(id[[i]], mapped_seqs)

            entrez_id[[i]] <- na.omit(unique(as.list(unlist( ref2eg[id[[i]]] ) )))
            nonunique_id[[i]] <- length(id[[i]]) - length(entrez_id[[i]])
            writeLines(paste0("Number of non-unique genes in cluster", n, " = ", nonunique_id[[i]]))
        }
        writeLines("converted gene symbols to entrez ids")
    }
    
    #### Create parameters for each cluster and each GO class
    parameters <- list(NULL)
    result <- vector("list", length = length(cluster_id)) 
        # Each cluster contains 5 lists 
        # (1 KEGG + 1 KEGGHyperGParams object)
        # (3 GO classes + 4 HyperGOresults objects per class)

    for (i in 1:length(cluster_id)) {
        for (j in 1:length(go.classes)) {
            writeLines(go.classes[j])
            n <- cluster_id[i]
            if (length(entrez_id[[i]]) > 0) {
                
                ### KEGG ###
                if (go.classes[j] == "KEGG") {
                    parameters <- new("KEGGHyperGParams",
                        geneIds=entrez_id[[i]],
                        universeGeneIds=universe, # can subset further
                        annotation=hs.annotations,
                        categoryName=go.classes[j],
                        pvalueCutoff=1, # no cutoff 
                        testDirection="over")
                } else {
                ### GO ###
                    parameters <- new("GOHyperGParams",
                        geneIds=entrez_id[[i]],
                        universeGeneIds=universe, # can subset further
                        annotation=hs.annotations,
                        ontology=go.classes[j],
                        pvalueCutoff=1, # no cutoff 
                        conditional=F, 
                        testDirection="over")
                }
                #### Hypergeometric tests ####
                if (exists("output")) { rm(output) } #endif
                output <- hyperGTest(parameters)
                result[[i]][[no.classes+1]][[j]] <- output # save object
                output <- summary(output)

                #### Multiple-Test Correction #### 
                if(nrow(output) == 0){
                    # no correction
                    result[[i]][[j]] <- NULL
                } else {
                    # Adjust p-value in new column
                    output$p.adjusted <- p.adjust(output[,2], "BH") #"fdr"
                    #re-arrange p.adjusted column next to pvalue column
                    output <- output[1:5,c(1:2, ncol(output), 3:(ncol(output)-1))] 
                    result[[i]][[j]]<- output
                }
            } else {
                result[[i]][[j]] <- NULL
                writeLines(
                    paste0("No differentially expressed genes in cluster ",
                           n,". no GO computed."))
            }
        }
        writeLines(
            paste0("Cluster: ", n,
                   " gene ontology analysis complete."))
    }
    
    names(result) <- cluster_id
    return(result)

}

#### Create GO Table : Each cluster has it's own column indicating the -log10 adjusted p-values

makeGOtable <- function(GOresults) {
    
    ## create dataframe of go terms
    terms <- Term(GOTERM)
    terms.df <- data.frame(ID = names(terms), GOname = terms, row.names = NULL)
    
    ## store metadata of GOresults
    numclusters <- length(GOresults)
    cat("Number of clusters: ", numclusters, "\n")
    cluster_names <- names(GOresults)
    cat("Cluster names: ", cluster_names, "\n")
    
    ## create dataframe to store results and temporary vector variable
    GOtable <- data.frame(ID = character())
    x <- vector()
    
    ## perform multiple-test correction on all 3 GO categories and add to results data.frame
    start.time <- Sys.time()
    for (cluster_id in 1:numclusters) {
        x <- pvalues(GOresults[[ cluster_names[cluster_id] ]][[4]][[1]])
        x <- append(x, pvalues(GOresults[[cluster_names[cluster_id]]][[4]][[2]]), after = length(x))
        x <- append(x, pvalues(GOresults[[cluster_names[cluster_id]]][[4]][[3]]), after = length(x))
        cat("Number of GO IDs in cluster ", cluster_names[cluster_id], " : ", length(x), "\n")
        x <- -log10(p.adjust(x, "BH")) 

        GOtable
        cat("Completed multiple test correction", "\n")
        df <- setNames(data.frame(names(x), x), c('ID', cluster_names[cluster_id]))

        GOtable <- merge(GOtable, df, by = 'ID', all=TRUE)
        cat("GO results for cluster", cluster_names[cluster_id], "complete\n")
    }
    end.time <- Sys.time()
    cat("Time elapsed in for loop: ", end.time - start.time, "\n")
    
    # replace NA with 0  
    GOtable[is.na(GOtable)] <- 0
    
    # left merge GO table with terms.df (all GO terms), results of GOtable are kept even if there are no matching descriptions in terms.df
    GOtable <- merge(GOtable, terms.df, by = 'ID', all.x = TRUE)
    
    # move adjusted GO description next to the GO.ID.
    GOtable <- GOtable[,c(1,length(GOtable), 2:(length(GOtable)-1))]
    
    return(GOtable)
}

# Create csv file summarizing all GO terms from all clusters, including columns for cluster, GO class, GO IDs, and GO terms, adjusted p-values and -log10(p.adjust), and Number of genes per GO term
makeKEGGsummary <- function(GOresults) {
    result <- data.frame(KEGGID = character(), Term = character())
    df <- NA
    kegg.locn <- 4
    
    ## store metadata of GOresults
    numclusters <- length(GOresults)
    cat("Number of clusters: ", numclusters, "\n")
    cluster_names <- names(GOresults)
    cat("Cluster names: ", cluster_names, "\n")
    
    for (cluster_id in 1:numclusters) {
        cat("Beginning summary of cluster ", cluster_names[cluster_id] ,"\n")
        df <- summary(GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[kegg.locn]])
        df$cluster <- cluster_names[cluster_id]
        df$Class <- GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[kegg.locn]]@testName
        df$p.adj <- p.adjust(df[,2],"BH") # FDR correction
        df$nlog10p.adj <- -log10(df$p.adj)
        df <- df[,c(1, (ncol(df)-4), (ncol(df)-2), 2, (ncol(df)-1):ncol(df), 3:6, (ncol(df)-3) )]
        colnames(df) <- c("KEGGID", colnames(df)[2:length(df)])
        if (cluster_id == 1) {
            result <- df
        } else {
            result <- rbind(result, df)
        }
        cat("KEGG computed","\n")
    }
    result[is.na(result$nlog10p.adj)] <- 0     # In column -log10(p.adjust) replace NA with 0  
    return(result)
}

makeGOsummary <- function(GOresults) {
    
    result <- data.frame(GOID = character(), Term = character())
    df <- NA
    
    ## store metadata of GOresults
    numclusters <- length(GOresults)
    cat("Number of clusters: ", numclusters, "\n")
    cluster_names <- names(GOresults)
    cat("Cluster names: ", cluster_names, "\n")
    
    for (cluster_id in 1:numclusters) {
        if (cluster_id == 1) {
            cat("Beginning summary of cluster ", cluster_names[cluster_id] ,"\n")
            df <- summary(GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[1]])
            df$cluster <- cluster_names[cluster_id]
            df$Class <- GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[1]]@testName[2]
            df$p.adj <- p.adjust(df[,2],"BH") # FDR correction
            df$nlog10p.adj <- -log10(df$p.adj)

            df <- df[,c(1, (ncol(df)-4), (ncol(df)-2), 2, (ncol(df)-1):ncol(df), 3:6, (ncol(df)-3) )]
            colnames(df) <- c("GOID", colnames(df)[2:length(df)])
            result <- df
            cat("BP computed","\n")
        } else {
            cat("Beginning summary of cluster ", cluster_names[cluster_id] ,"\n")
            df <- summary(GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[1]])
            df$cluster <- cluster_names[cluster_id]
            df$Class <- GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[1]]@testName[2]
            df$p.adj <- p.adjust(df[,2],"BH") # FDR correction
            df$nlog10p.adj <- -log10(df$p.adj)
            
            df <- df[,c(1, (ncol(df)-4), (ncol(df)-2), 2, (ncol(df)-1):ncol(df), 3:6, (ncol(df)-3) )]
            colnames(df) <- c("GOID", colnames(df)[2:length(df)])
            result <- rbind(result, df)
            cat("BP computed","\n")
        }
        
        df <- summary(GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[2]])
        df$cluster <- cluster_names[cluster_id]
        df$Class <- GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[2]]@testName[2]
        df$p.adj <- p.adjust(df[,2],"BH") # FDR correction
        df$nlog10p.adj <- -log10(df$p.adj)

        df <- df[,c(1, (ncol(df)-4), (ncol(df)-2), 2, (ncol(df)-1):ncol(df), 3:6, (ncol(df)-3) )]
        colnames(df) <- c("GOID", colnames(df)[2:length(df)])
        result <- rbind(result, df)
        cat("MF computed","\n")

        df <- summary(GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[3]])
        df$cluster <- cluster_names[cluster_id]
        df$Class <- GOresults[[ cluster_names[cluster_id] ]][[no.classes+1]][[3]]@testName[2]
        df$p.adj <- p.adjust(df[,2],"BH") # FDR correction
        df$nlog10p.adj <- -log10(df$p.adj)

        df <- df[,c(1, (ncol(df)-4), (ncol(df)-2), 2, (ncol(df)-1):ncol(df), 3:6, (ncol(df)-3) )]
        colnames(df) <- c("GOID", colnames(df)[2:length(df)])
        result <- rbind(result, df)
        cat("CC computed","\n") 

    }
    
    result[is.na(result$nlog10p.adj)] <- 0     # In column -log10(p.adjust) replace NA with 0  
    return(result)
}

### Function to create GO bar plots
GOBarPlots <- function(GOsum, p.threshold = 0.05, topn = 10){
    
    cat("Threshold for p.adj: ", p.threshold, 
        "\nNumber of Top GO terms:", topn, "\n")
    
    ## store metadata of GOsum
    numclusters <- length(unique(GOsum$cluster))
    cat("Number of clusters: ", numclusters, "\n")
    cluster_names <- as.character(unique(GOsum$cluster))
    cat("Cluster names: ", cluster_names, "\n")
    
    ## Initialize output object
    result <- vector(mode = "list", length = numclusters)
    names(result) <- cluster_names
    
    ## Create plots
    for (cluster_id in 1:numclusters) {
        n <- cluster_names[cluster_id]
        # cat("Plotting for Cluster: ", n, "\n")
        
        ## subset GOsum to reduce GOIDs for the specific cluster
        temp <- subset(GOsum, (p.adj < p.threshold & cluster == n))
        
        ## Select top_n (10) Go terms.
        temp_split <- split(temp, temp$Class)
        # writeLines(length(temp_split))
        
        if (!is.null(temp_split$BP)) {
            top_BP <- top_n(temp_split$BP, topn, wt=nlog10p.adj)
        } else {
            top_BP <- data.frame(matrix(ncol = ncol(temp), nrow = 0))
            colnames(top_BP) <- colnames(temp)
        }
        if (!is.null(temp_split$MF)) {
            top_MF <- top_n(temp_split$MF, topn, wt=nlog10p.adj)
        } else {
            top_MF <- data.frame(matrix(ncol = ncol(temp), nrow = 0))
            colnames(top_MF) <- colnames(temp)
        }
        if (!is.null(temp_split$CC)) {
            top_CC <- top_n(temp_split$CC, topn, wt=nlog10p.adj)
        } else {
            top_CC <- data.frame(matrix(ncol = ncol(temp), nrow = 0))
            colnames(top_CC) <- colnames(temp)
        }

        temp <- rbind(top_BP, top_MF, top_CC)
        temp$GOID <- factor(temp$GOID, levels=rev(temp$GOID))
        temp$Term <- factor(temp$Term, levels=rev(temp$Term))
        temp$Class <- factor(temp$Class, levels=unique(temp$Class))            
        
        ## Create the plot
        p <- ggplot(temp, aes(x=nlog10p.adj, y=Term, fill=Class)) + 
                geom_bar(stat='identity') + 
                scale_fill_manual("GO Class", values = c("BP" = "slateblue3", "MF" = "lightseagreen", "CC" = "darkgoldenrod1")) +
                geom_text(aes(label=Count), color="azure4", hjust = -0.5) +
                ggtitle(paste0("Cluster ",n)) + theme(plot.title = element_text(hjust = 0.5))
        
        ## save the plot in the output object
        result[[n]] <- p
    }
    cat("Plots created")
    return(result)
}

### Function to return lists of Gene symbols for all GO terms in a cluster
### Input: List of GOHyperGResult objects (BP, MF, CC)
### Output: list of GO terms and the associated gene symbols
genesbyGO <- function(GOResult, cluster){ 
    eg2ref <- org.Hs.egSYMBOL
    result <- data.frame(ID = character(), 
                         # Term = character(), 
                         Genes = character(),
                         Class = character(),
                         cluster = character(),
                         row.names = NULL)
    for (j in 1:length(go.classes)) { # Loop through each Ontology class
        writeLines(paste0("Processing Ontology Category: ", go.classes[j]))
        genesbycat <- geneIdsByCategory(GOResult[[j]])
        writeLines(paste0("Number of GO terms in this Ontology Class: ",length(genesbycat)))
        
        for (i in 1:length(genesbycat)){ # Loop through each GO/KEGG
            geneSymbols <- na.omit(unique(as.list(unlist( eg2ref[genesbycat[[i]]] ) ) ))
            temp <- do.call(cbind, list(names(genesbycat)[i], 
                                        geneSymbols,
                                        GOResult[[j]]@testName[2],
                                        cluster-1
                                        ))
            if (j == 1 & i == 1) {
                result <- temp
                colnames(result) <- c("ID", "Genes", "Class", "cluster")
            } else {
                result <- rbind(result, temp)
            }
        }
    }
    writeLines(paste0("Genes Extracted for cluster: ", cluster ))
    return(result)
}

### Convert entrez ids of all genes back to gene symbols using lapply loop on the GO classes for each cluster
### Input: list resulting from GO analysis (nested list of clusters and their summary dataframes and GOHyperGResult object)
### Output: nested list of clusters and their Ontology classes, GO terms and associated gene symbols. 
allGOgenes <- function (GOresults) {
    results <- NULL
    for (i in 1:length(GOresults)) {
        writeLines(paste0("Processing cluster ", i))
        # class(GOresults[[i]][[no.classes+1]])
        if (i == 1) {
            results <- genesbyGO(GOresults[[i]][[no.classes+1]], i)
        } else {
            results <- rbind(results, 
                            genesbyGO(GOresults[[i]][[no.classes+1]], i))
        }
    }
    writeLines("Genes Extracted for all clusters")
    return(results)
}
