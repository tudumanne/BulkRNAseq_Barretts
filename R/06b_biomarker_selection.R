#This script contains the custom R code used for biomarker filtering based on following criteria
#log2FoldChange>2, p.adj<0.05 and median TPM>5
#Input files: summarized count data and DEG lists from both simple and adjusted models

#load required packages
library(dplyr)

#pair-wise comparisons used for biomarker selection
comparisons <- c(
  "LGD_vs_NDBE",
  "HGD_vs_LGD",
  "EAC_vs_HGD",
  "HGD_vs_NDBE",
  "EAC_vs_NDBE",
  "all_vs_NDBE", #LGD, HGD and EAC
  "dys_vs_NDBE" #LGD and HGD
)

#gene annotations
genes <- read.csv("tx2gene.gencode.v36.pa_edited.csv", header = TRUE)

#load summarized count data
load("counts_summary.RData")
load("tpm_summary.RData")

#X and Y used as generalized names below in relevance to running the 'for' loop
for (comp in comparisons) {
  
  #extract X and Y from pair-wise comparison name
  X <- strsplit(comp, "_vs_")[[1]][1]
  Y <- strsplit(comp, "_vs_")[[1]][2]
  
  #read DEG list - simple model
  deg_simple <- read.csv(
    paste0(comp, "_sigFC1.csv"),
    header = TRUE
  )
  #read DEG list - adjusted model
  deg_adjusted <- read.csv(
    paste0("adjusted_", comp, "_sigFC1.csv"),
    header = TRUE
  )
  
  #rename columns in adjusted model DEG lists
  deg_adjusted <- deg_adjusted %>% 
    rename(
      log2FoldChange.adjusted = log2FoldChange,
      padj.adjusted = padj
    )
  
  #remove unecessary columns
  deg_simple2 <- deg_simple[, !names(deg_simple) %in% 
                              c("baseMean", "lfcSE", "stat", "pvalue")]
  
  deg_adjusted2 <- deg_adjusted[, !names(deg_adjusted) %in% 
                                  c("baseMean", "lfcSE", "stat", "pvalue")]
  
  #overlap of simple and adjusted model DEG lists
  deg_overlap <- deg_simple2 %>% 
    inner_join(deg_adjusted2, by = "X") %>%
    rename(gene_id = X)
  
  #add gene annotation
  deg_overlap2 <- deg_overlap %>% 
    inner_join(genes, by = "gene_id")
  
  #extract average normalized counts from DESeq2
  X_avg_counts <- paste0(X, "_avg_counts")
  Y_avg_counts <- paste0(Y, "_avg_counts")
  
  #merge with overlap list
  deg_overlap2 <- deg_overlap2 %>% 
    left_join(
      select(
        counts_summary,
        gene_name,
        all_of(Y_avg_counts),
        all_of(X_avg_counts)
      ),
      by = "gene_name"
    )
  
  #rename above selected columns to standard names
  deg_overlap2 <- deg_overlap2 %>%
    rename(
      X_avg_counts = all_of(X_avg_counts),
      Y_avg_counts = all_of(Y_avg_counts)
    )
  
  #calculate ratio
  deg_overlap2$XY_ratio_counts <- 
    deg_overlap2$X_avg_counts / deg_overlap2$Y_avg_counts

  #extract median TPM columns
  X_median_tpm <- paste0(X, "_median_tpm")
  Y_median_tpm <- paste0(Y, "_median_tpm")
  
  #merge
  deg_overlap2 <- deg_overlap2 %>% 
    left_join(
      select(
        tpm_summary,
        gene_id,
        all_of(Y_median_tpm),
        all_of(X_median_tpm)
      ),
      by = "gene_id"
    )
  
  #rename to standard names
  deg_overlap2 <- deg_overlap2 %>%
    rename(
      X_median_tpm = all_of(X_median_tpm),
      Y_median_tpm = all_of(Y_median_tpm)
    )
  
  #final combined table 
  XY_deg_overlap <- select(
    deg_overlap2,
    gene_id,
    gene_name,
    log2FoldChange,
    padj,
    log2FoldChange.adjusted,
    padj.adjusted,
    Y_avg_counts,
    X_avg_counts,
    XY_ratio_counts,
    Y_median_tpm,
    X_median_tpm
  )
  
  #define and extract up-regulated gene shortlist
  XY_subset_up_shortlist <- XY_deg_overlap %>% 
    filter(
      log2FoldChange >= 2,
      XY_ratio_counts >= 5,
      X_median_tpm >= 5
    )
  
  #define and extract down-regulated gene shortlist
  XY_subset_down_shortlist <- XY_deg_overlap %>% 
    filter(
      log2FoldChange <= -2,
      XY_ratio_counts <= 0.2,
      X_median_tpm == 0
    )
  
  #export as .csv files
  write.csv(
    XY_subset_up_shortlist,
    paste0(comp, "_upregulated_shortlist.csv"),
    row.names = FALSE
  )
  
  write.csv(
    XY_subset_down_shortlist,
    paste0(comp, "_downregulated_shortlist.csv"),
    row.names = FALSE
  )
  
  #save .RData file
  save(
    XY_deg_overlap,
    XY_subset_up_shortlist,
    XY_subset_down_shortlist,
    file = paste0(comp, ".RData")
  )
  
}