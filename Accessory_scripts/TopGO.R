# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("topGO")
library(topGO)

############## ############## ############## ############## ############## ############## ############## ############## 
############## extracting expression data for GO terms

genesInTerm(Tpa_GODATA_BP, "GO:0051321")

### Function to get transcript IDs annotated with particular GO terms

get_trans_from_GOs = function(GOLIST,GODATA){
  #mygenes <- genesInTerm(GODATA, GOLIST)
  out_me <- list()
  
  for (i in 1:length(GOLIST))
  {
    print(GOLIST[i])	
    mygenes <- genesInTerm(GODATA, GOLIST[i])
    #print(mygenes)
    mygenes_GOs <- append(GOLIST[i], mygenes )
    print(mygenes_GOs)
    out_me = append(out_me, mygenes)
  }
  
  print(out_me)
  
  all_genes <- c()
  for (n in seq(1, length(out_me))){
    print(n)
    all_genes <- c(all_genes, out_me[[n]])
  }	
  
  all_genes_u <- unique(all_genes)  
  return(all_genes_u)
  
}


######################################
######### load annotations


setwd("/Users/drp22jhz/Documents/University/Lausanne/Timema_transcriptomes/Timema_DC_MSCI_code/data/annot") 

## Droso annotated
geneID2GO_Tpa_DROSO <- readMappings(file = "Tpa_LRv5a_mtDNAv350_v2.2_droso_b2g.gff_droso_ontology_term_fortopgo.txt")


## take a look
str(head(geneID2GO_Tpa_DROSO))
head(geneID2GO_Tpa_DROSO)


geneID2GO_Tpa_DROSO$Tpa_000270
geneID2GO_Tpa_DROSO$Tpa_028282


### then need to get the genes as a named numeric vector
RGL_Tpa <- as.data.frame(read.table("Tpa_LRv5a_mtDNAv350_v2.2_droso_b2g.gff_droso_ontology_term_fortopgo_genelist.txt"))
RGL_Tpa$V2 <- rep(1, length(RGL_Tpa[,1])) ### not doing enrichment - so just 1s
RGL_Tpa_1_GL <- RGL_Tpa$V2
names(RGL_Tpa_1_GL) <- RGL_Tpa$V1

subset(RGL_Tpa, RGL_Tpa$V1 == "Tpa_000270")


#### from here you can make the GO hierarchy

### make rule for classing sig / non-sig (for fisher test)
topDiffGenes <- function(allScore) {return(allScore < 0.05)}

#### make GOdata object
#### specifing like this make the GO universe from the genes in the DJP_geneList 
#### setting node size as 5 so at least 5 genes must be annot per GO terms 
#### do enrichment tests

Tpa_GODATA_BP = new("topGOdata", ontology = "BP", allGenes = RGL_Tpa_1_GL, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = geneID2GO_Tpa_DROSO, nodeSize = 1) ### note just BP here (also can use MF and CC)


############## ############## ############## ############## ############## ############## ############## ############## 
############## extracting expression data for GO terms

### get genes with immune terms annotated
meiosis_cell_cycle_terms = c("GO:0051321") 


## This will give the geneIDs of the genes annotated with the GO terms you specified
Tpa_meiosis_cell_cycle_genes <- get_trans_from_GOs(meiosis_cell_cycle_terms, Tpa_GODATA_BP) ## will give warnings for GOs not found



############## ############## ############## ############## ############## ############## ############## ############## 
############## 
############## ############## ############## ############## ############## ############## ############## ############## 
############## 
## Droso annotated
geneID2GO_Tps_DROSO <- readMappings(file = "Tps_LRv5b_mtDNAv350_v2.2_droso_b2g.gff_droso_ontology_term_fortopgo.txt")


## take a look
str(head(geneID2GO_Tps_DROSO))
head(geneID2GO_Tps_DROSO)


geneID2GO_Tps_DROSO$Tps_000270
geneID2GO_Tps_DROSO$Tps_028282


### then need to get the genes as a named numeric vector
RGL_Tps <- as.data.frame(read.table("Tps_LRv5b_mtDNAv350_v2.2_droso_b2g.gff_droso_ontology_term_fortopgo_genelist.txt"))
RGL_Tps$V2 <- rep(1, length(RGL_Tps[,1])) ### not doing enrichment - so just 1s
RGL_Tps_1_GL <- RGL_Tps$V2
names(RGL_Tps_1_GL) <- RGL_Tps$V1

subset(RGL_Tps, RGL_Tps$V1 == "Tps_000270")


#### from here you can make the GO hierarchy

### make rule for classing sig / non-sig (for fisher test)
topDiffGenes <- function(allScore) {return(allScore < 0.05)}

#### make GOdata object
#### specifing like this make the GO universe from the genes in the DJP_geneList 
#### setting node size as 5 so at least 5 genes must be annot per GO terms 
#### do enrichment tests

Tps_GODATA_BP = new("topGOdata", ontology = "BP", allGenes = RGL_Tps_1_GL, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = geneID2GO_Tps_DROSO, nodeSize = 1) ### note just BP here (also can use MF and CC)


############## ############## ############## ############## ############## ############## ############## ############## 
############## extracting expression data for GO terms

### get genes with immune terms annotated
meiosis_cell_cycle_terms = c("GO:0051321") 


## This will give the geneIDs of the genes annotated with the GO terms you specified
Tps_meiosis_cell_cycle_genes <- get_trans_from_GOs(meiosis_cell_cycle_terms, Tps_GODATA_BP) ## will give warnings for GOs not found





############## ############## ############## ############## ############## ############## ############## ############## 
##############
############## ############## ############## ############## ############## ############## ############## ############## 
##############
############## ############## ############## ############## ############## ############## ############## ############## 
##############
## Droso annotated
geneID2GO_Tce_DROSO <- readMappings(file = "Tce_LRv5a_mtDNAv350_v2.2_droso_b2g.gff_droso_ontology_term_fortopgo.txt")


## take a look
str(head(geneID2GO_Tce_DROSO))
head(geneID2GO_Tce_DROSO)


geneID2GO_Tce_DROSO$Tce_000270
geneID2GO_Tce_DROSO$Tce_028282


### then need to get the genes as a named numeric vector
RGL_Tce <- as.data.frame(read.table("Tce_LRv5a_mtDNAv350_v2.2_droso_b2g.gff_droso_ontology_term_fortopgo_genelist.txt"))
RGL_Tce$V2 <- rep(1, length(RGL_Tce[,1])) ### not doing enrichment - so just 1s
RGL_Tce_1_GL <- RGL_Tce$V2
names(RGL_Tce_1_GL) <- RGL_Tce$V1

subset(RGL_Tce, RGL_Tce$V1 == "Tce_000270")


#### from here you can make the GO hierarchy

### make rule for classing sig / non-sig (for fisher test)
topDiffGenes <- function(allScore) {return(allScore < 0.05)}

#### make GOdata object
#### specifing like this make the GO universe from the genes in the DJP_geneList 
#### setting node size as 5 so at least 5 genes must be annot per GO terms 
#### do enrichment tests

Tce_GODATA_BP = new("topGOdata", ontology = "BP", allGenes = RGL_Tce_1_GL, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = geneID2GO_Tce_DROSO, nodeSize = 1) ### note just BP here (also can use MF and CC)


############## ############## ############## ############## ############## ############## ############## ############## 
############## extracting expression data for GO terms

### get genes with immune terms annotated
meiosis_cell_cycle_terms = c("GO:0051321") 


## This will give the geneIDs of the genes annotated with the GO terms you specified
Tce_meiosis_cell_cycle_genes <- get_trans_from_GOs(meiosis_cell_cycle_terms, Tce_GODATA_BP) ## will give warnings for GOs not found


##### output 

write.csv(c(Tce_meiosis_cell_cycle_genes, Tps_meiosis_cell_cycle_genes, Tpa_meiosis_cell_cycle_genes), "meiosis_cell_cycle_genes.csv", row.names = F)





