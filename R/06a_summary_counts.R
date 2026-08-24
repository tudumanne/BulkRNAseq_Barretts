#This script contains the custom R code used to summarize count data required for biomarker selection
#Input files: normalized count data from DESeq2 (simple model) and TPM counts from Salmon

# Load library
library("dplyr")

#normalized counts from DESeq2
counts = read.csv("normalized_counts.csv", header=TRUE)
#rename columns
counts <- counts %>% 
  rename("gene_id"="V2", "gene_name"="V3")

#subset counts as required here based on sample ID and calculate average
#B=NDBE, L=LGD, H=HGD, C=EAC

subset_counts_NDBE = counts %>% select(gene_id | starts_with("B"))
subset_counts_NDBE$NDBE_avg_counts = rowMeans(subset_counts_NDBE[c(2:30)])
counts_NDBE = subset_counts_NDBE %>%
  select(gene_id, NDBE_avg_counts)

subset_counts_LGD = counts %>% 
  select(gene_id | starts_with("L"))
subset_counts_LGD$LGD_avg_counts = rowMeans(subset_counts_LGD[c(2:22)])
counts_LGD = subset_counts_LGD %>%
  select(gene_id, LGD_avg_counts)

subset_counts_HGD = counts %>% 
  select(gene_id | starts_with("H"))
subset_counts_HGD$HGD_avg_counts = rowMeans(subset_counts_HGD[c(2:12)])
counts_HGD = subset_counts_HGD %>%
  select(gene_id, HGD_avg_counts)

subset_counts_EAC= counts %>% select(gene_id | starts_with("C"))
subset_counts_EAC$EAC_avg_counts = rowMeans(subset_counts_EAC[c(2:15)])
counts_EAC = subset_counts_EAC %>%
  select(gene_id, EAC_avg_counts)

#here 'all' refers to LGD, HGD and EAC
subset_counts_all = counts %>% select(gene_id | starts_with("L")| starts_with("H")| starts_with("C"))
subset_counts_all$all_avg_counts = rowMeans(subset_counts_all[c(2:47)])
counts_all = subset_counts_all %>%
  select(gene_id, all_avg_counts)

#here 'dys' refers to LGD and HGD
subset_counts_dys = counts %>% select(gene_id | starts_with("L")| starts_with("H"))
subset_counts_dys$dys_avg_counts = rowMeans(subset_counts_dys[c(2:33)])
counts_dys = subset_counts_dys %>%
  select(gene_id, dys_avg_counts)

counts1 = counts_NDBE %>% inner_join(counts_LGD, by = "gene_id")
counts2 = counts1 %>% inner_join(counts_HGD, by = "gene_id")
counts3 = counts2 %>% inner_join(counts_EAC, by = "gene_id")
counts4 = counts3 %>% inner_join(counts_all, by = "gene_id")
counts_summary = counts %>% inner_join(counts_dys, by = "gene_id")

write.csv(counts_summary, "counts_summary.csv")
save(counts_summary, file="counts_summary.RData")

salmon_tpm = read.csv("salmon_counts.csv", header=TRUE)

#subset counts as required here based on sample ID
#B=NDBE, L=LGD, H=HGD, C=EAC
salmon_tpm$NDBE_median_tpm = apply(salmon_tpm[,grep("B",colnames(salmon_tpm))], 1, median, na.rm=T)
salmon_tpm$LGD_median_tpm = apply(salmon_tpm[,grep("L",colnames(salmon_tpm))], 1, median, na.rm=T)
salmon_tpm$HGD_median_tpm = apply(salmon_tpm[,grep("H",colnames(salmon_tpm))], 1, median, na.rm=T)
salmon_tpm$EAC_median_tpm = apply(salmon_tpm[,grep("C",colnames(salmon_tpm))], 1, median, na.rm=T)
all = c("L", "H", "C")
salmon_tpm$all_median_tpm = apply(salmon_tpm[,grep(paste(all,collapse="|"),colnames(salmon_tpm))], 1, median, na.rm=T)
dys = c("L","H")
salmon_tpm$dys_median_tpm = apply(salmon_tpm[,grep(paste(dys,collapse="|"),colnames(salmon_tpm))], 1, median, na.rm=T)

tpm_summary = salmon_tpm %>%
  select(gene_id, NDBE_median_tpm, LGD_median_tpm, HGD_median_tpm, EAC_median_tpm,
         all_median_tpm, dys_median_tpm)

write.csv(tpm_summary, "tpm_summary.csv")
save(tpm_summary, file="tpm_summary.RData")
