#This script run the formatting required for decoupleR analyses by adding gene names to DEG lists from DESeq2 analysis

#load required libraries
library(dplyr)

#define DEG lists of interest
deg_lists <- list(
  "NDBE_vs_NSq_allgenes.csv",
  "LGD_vs_NDBE_allgenes.csv",
  "HGD_vs_LGD_allgenes.csv",
  "EAC_vs_HGD_allgenes.csv"
)

#gene annotations
genelist = read.csv("tx2gene.gencode.v36.pa_edited.csv", header=TRUE)

#read in DEG lists, merge gene names, save as .csv files
for (file in deg_lists){
  deg_list <-read.csv(file, header=TRUE)
  
  deg_list_genename <- deg_list %>% left_join(genelist, by = "gene_id")
  
  output_file <- sub(
    "_allgenes\\.csv$",
    "allgenes_genename.csv", 
    file
  )
  
  write.csv(
    deg_list_genename, 
    output_file
  )
}
