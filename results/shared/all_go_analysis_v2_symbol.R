all_go_analysis_v2_symbol <- function(gene, org,...){
	     library(GOstats)
	     library(KEGG.db)

	     class <- c("KEGG","BP","MF","CC")
	     if(org == "human"){
	     	   library(org.Hs.eg.db)
		   ann <- "org.Hs.eg.db"
	     	   ref2eg <- org.Hs.egSYMBOL2EG
		   mapped_seqs <- mappedkeys(ref2eg)
		   universe <- ls(org.Hs.egSYMBOL)
	     }else if(org == "mouse"){
	     	   library(org.Mm.eg.db)
		   ann <- "org.Mm.eg.db"
	     	   ref2eg <- org.Mm.egSYMBOL2EG
		   mapped_seqs <- mappedkeys(ref2eg)
		   universe <- ls(org.Mm.egSYMBOL)
	     }else{
		exit
	     }
	     
	     #get cluster number
	     id <- unique(as.character(gene))
	     #print(id)
	     if(length(id) == 0){
	     		   return(NULL)
	     }

	     id <- intersect(id, mapped_seqs)
	     entrezid <- na.omit(unique(unlist(as.list(ref2eg[id]))))

		   #GO
		   if(org=="human"){
			 GOid <- mget(entrezid, org.Hs.egGO, ifnotfound=NA)
		   }else if(org == "mouse"){
		   	 GOid <- mget(entrezid, org.Mm.egGO, ifnotfound=NA)
		   }

		   result <- vector("list",length(class))
		   for(j in 1:length(class)){
		   print(class[j])
		   #KEGG
		   if(class[j] == "KEGG"){
		   para <- new("KEGGHyperGParams",
		   geneIds=entrezid,
		   universeGeneIds=universe,
		   annotation=ann,
		   categoryName=class[j],
		   pvalueCutoff=1,
		   testDirection="over")
		   }else{
		   haveid <- sapply(GOid, function(x){
                         if(length(x) == 1 && is.na(x)){
                                      FALSE
                         }else{
                                      onts <- Biobase::subListExtract(x, "Ontology",simplify=FALSE)
                                      class[j] %in% onts
                         }})
		   if(length(haveid[haveid]) == 0){
		   			     result[[j]] <- NULL
		   			     next
		   }
		   entrezid <- names(haveid[haveid])
		   #print(entrezid)

		   para <- new("GOHyperGParams",
		   geneIds=entrezid,
		   universeGeneIds=universe,
		   annotation=ann,
		   ontology=class[j],
		   pvalueCutoff=1,
		   conditional=F,
		   testDirection="over")
		   }
		   data <- hyperGTest(para)
		   data <- summary(data)
		   
		   data <- data[data[,6]<500,]
		   data <- data[data[,6]>50,]
		   #data <- data

		   if(nrow(data) == 0){
		   result[[j]] <- NULL
		   }else{
		   data[,2] <- p.adjust(data[,2],"fdr")
		   result[[j]] <- data
		   #result[[j]] <- multi_correct(data, "BH", 0.05)
		   }
		   }

	     return(result)
}