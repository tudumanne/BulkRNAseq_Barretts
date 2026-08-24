#Supplementary figure S2C
#This script contains the R code used to run ORA analysis of GO terms (Biological Processes)
#Input files: 
#DEG lists containing significant results (all disease stages vs NSq, |log2FC|>1 and adj.p<0.05) 
#and a list of all genes as the universe

#References:
#https://yulab-smu.top/biomedical-knowledge-mining-book/

#load required packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggnewscale)
library(DOSE)
library(ReactomePA)

#read in gene list 
data1 = read.csv("all_vs_NSq_up_sig.csv", header = TRUE)
data2 = read.csv("all_vs_NSq_down_sig.csv", header = TRUE)

#set the desired organism
organism = org.Hs.eg.db

#gene list used as the universe
genelist = read.csv("all_genes.csv", header = TRUE)

gene_list = genelist$genes

gene_list1 = data1$gene_id
gene_list2 = data2$gene_id

geneClusters <- list(
  all_up = gene_list1,
  all_down = gene_list2
)

#run ORA
cc_ora <- compareCluster(
  geneCluster  = geneClusters,
  fun          = "enrichGO",
  universe=gene_list,
  ont ="BP", 
  keyType = "ENSEMBL", 
  pvalueCutoff = 0.05,
  OrgDb = organism, 
  pAdjustMethod = "none"
)

head(as.data.frame(cc_ora))
#dotplot
dotplot(cc_ora, showCategory = 10) #+ facet_grid(.~ Cluster)
