#Figure 2A
#This script contains the R code used to run GSEA analysis of GO terms (Biological Processes)
#Input files: DEG lists containing all genes

#References:
#https://yulab-smu.top/biomedical-knowledge-mining-book/

#load required packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggnewscale)
library(DOSE)
library(ReactomePA)

file_names <- c("NDBE_vs_NSq_allgenes.csv", "LGD_vs_NDBE_allgenes.csv", "HGD_vs_LGD_allgenes.csv", "EAC_vs_HGD_allgenes.csv")
list_names  <- c("nsq_ndbe_list", "lgd_ndbe_list", "hgd_lgd_list", "eac_hgd_list")

for (i in 1:length(file_names)) {
  
  #read the DEG lists
  temp_data <- read.csv(file_names[i], header=TRUE)
  
  #extract log2FC and assign names from column X
  temp_vector <- temp_data$log2FoldChange
  names(temp_vector) <- temp_data$X
  
  #clean NA values and sort by decreasing log2FC
  temp_vector <- na.omit(temp_vector)
  temp_vector <- sort(temp_vector, decreasing=TRUE)
  
  #save the lists
  assign(list_names[i], temp_vector)
}

#define gene clusters
geneClusters <- list(
  NDBE_NSq = nsq_ndbe_list,
  LGD_NDBE = lgd_ndbe_list,
  HGD_LGD = hgd_lgd_list,
  EAC_HGD = eac_hgd_list
)

#set the desired organism
organism = org.Hs.eg.db

#run GSEA 
cc_gsea <- compareCluster(
  geneCluster  = geneClusters,
  fun          = "gseGO",
  OrgDb        = organism,
  keyType      = "ENSEMBL",
  ont          = "BP",
  minGSSize = 10, 
  maxGSSize = 500, 
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  eps = 0
)

head(as.data.frame(cc_gsea))
#dotplot
dotplot(cc_gsea, showCategory = 6) + facet_grid(.~.sign)
