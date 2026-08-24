#This script contains the R code used to run differential gene expression analysis using DESeq2 (adjusted model)
#Input files: gene-level read count data from STAR, sample metadata including immune and stromal scores, and gene annotations

#References:
#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#load required packages
library(DESeq2)
library(ggplot2)
library(PCAtools)
library(ggrepel)

#sample metadata
sampletable <- read.table("sample_list_shortid.txt", header=T, sep="\t")

#set row names to match sample names
rownames(sampletable) <- sampletable$SampleName
head(sampletable)

#check the number of rows and the number of columns
nrow(sampletable)
ncol(sampletable)

#center and scale co-variates
sampletable$stromal.ESTIMATE = scale(sampletable$stromal.ESTIMATE, center=TRUE)
sampletable$immune.ESTIMATE = scale(sampletable$immune.ESTIMATE, center=TRUE)

#define the disease stages i.e. Classification as factors
Classification <- factor(c(rep("Nsq", 10), rep("NDBE", 29), rep("LGD", 21), rep("HGD", 11), rep("EAC", 14)))

#load count data from STAR, individual file per sample
#directory is the path to the folder where count data are stored
#design refers to the variables of interest, i.e. Classification/disease stage and co-variates
sampletable$Classification <- factor(sampletable$Classification, levels = c("NSq", "NDBE", "LGD", "HGD", "EAC"))
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampletable,
                                      directory = "raw_counts_all",
                                      design = ~ Classification + stromal.ESTIMATE + immune.ESTIMATE)
dds

#load annotation files
gene.annotation <- read.csv("tx2gene.gencode.v36.pa_edited.csv", header=TRUE)
#reorder and check if the annotation (gene ID) matches row names of dds object
gene.annotation <- gene.annotation[match(rownames(dds), gene.annotation$gene.id),]
all(rownames(dds) == gene.annotation$gene.id)
#check for NA values and duplicates in gene annotation
sum(is.na(gene.annotation$gene.id))
sum(duplicated(gene.annotation$gene.id))

#write raw counts to a .csv file
res <- as.data.frame(counts(dds, normalized=FALSE))
write.csv(res, "adjusted_raw_counts.csv", sep=",")

#number of genes in the dds object
nrow(dds)

#filter based on counts if required
#smallestGroupSize <- 10
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]
#number of genes left after low-count filtering:
#nrow(dds)

#run DESeq
design(dds) <- ~ 1 + Classification + stromal.ESTIMATE + immune.ESTIMATE
dds <- DESeq(dds)

#check the coefficients estimated by DEseq
resultsNames(dds)

#compute normalized counts (log2 transformed, + 1 to avoid errors)
norm_counts_untransformed <- counts(dds, normalized = TRUE)
norm_counts_log <- log2(counts(dds, normalized = TRUE)+1)

#add the gene symbols
norm_counts_untransformed_symbols <- merge(unique(gene.annotation[,2:3]), data.frame(ID=rownames(norm_counts_untransformed), norm_counts_untransformed), by=1, all=F)
norm_counts_log_symbols <- merge(unique(gene.annotation[,2:3]), data.frame(ID=rownames(norm_counts_log), norm_counts_log), by=1, all=F)
norm_counts_untransformed_symbols$V2 <- sub('\\.[0-9]*$', '', norm_counts_untransformed_symbols$V2)
norm_counts_log_symbols$V2 <- sub('\\.[0-9]*$', '', norm_counts_log_symbols$V2)

#write normalized counts to text files
write.table(norm_counts_untransformed_symbols, "adjusted_normalized_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")
write.table(norm_counts_log_symbols, "adjusted_normalized_counts_log.txt", quote=F, col.names=T, row.names=F, sep="\t")

#define the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat

#calculate coefficient vectors for each group
nsq <- colMeans(mod_mat[dds$Classification == "NSq", ])
ndbe <- colMeans(mod_mat[dds$Classification == "NDBE", ])
lgd <- colMeans(mod_mat[dds$Classification == "LGD", ])
hgd <- colMeans(mod_mat[dds$Classification == "HGD", ])
eac <- colMeans(mod_mat[dds$Classification == "EAC", ])

#obtain results for each pairwise contrast
res_ndbe_nsq <- results(dds, contrast = ndbe - nsq)
res_lgd_ndbe <- results(dds, contrast = lgd - ndbe)
res_hgd_lgd <- results(dds, contrast = hgd - lgd)
res_eac_hgd <- results(dds, contrast = eac - hgd)

res_lgd_nsq <- results(dds, contrast = lgd - nsq)
res_hgd_nsq <- results(dds, contrast = hgd - nsq)
res_eac_nsq <- results(dds, contrast = eac - nsq)

res_hgd_ndbe <- results(dds, contrast = hgd - ndbe)
res_eac_ndbe <- results(dds, contrast = eac - ndbe)

#define vector of coefficients groups of interest
ndbe_lgd_hgd_eac <- colMeans(mod_mat[dds$Classification %in% c("NDBE", "LGD", "HGD","EAC"),])
lgd_hgd <- colMeans(mod_mat[dds$Classification %in% c("LGD", "HGD"),])
lgd_hgd_eac <- colMeans(mod_mat[dds$Classification %in% c("LGD", "HGD","EAC"),])

#contrasts of interest
res_all_nsq <- results(dds, contrast = ndbe_lgd_hgd_eac - nsq)
res_all_ndbe <- results(dds, contrast = lgd_hgd_eac - ndbe)
res_dys_ndbe <- results(dds, contrast = lgd_hgd - ndbe)
res_eac_lgdhgd <- results(dds, contrast = eac - lgd_hgd)

#list of comparisons to export
deg_lists <- list(NDBE_vs_NSq = res_ndbe_nsq,
                  LGD_vs_NDBE = res_lgd_ndbe, 
                  HGD_vs_LGD = res_hgd_lgd, 
                  EAC_vs_HGD = res_eac_hgd, 
                  LGD_vs_NSq = res_lgd_nsq, 
                  HGD_vs_NSq = res_hgd_nsq, 
                  EAC_vs_NSq = res_eac_nsq, 
                  HGD_vs_NDBE = res_hgd_ndbe, 
                  EAC_vs_HGD = res_eac_ndbe, 
                  all_vs_NSq = res_all_nsq, 
                  all_vs_NDBE = res_all_ndbe, 
                  dys_vs_NDBE = res_dys_ndbe, 
                  EAC_vs_dys = res_eac_lgdhgd)

#save deg lists as .csv files
for(name in names(deg_lists)){
  write.csv(
    as.data.frame(deg_lists[[name]]),
    file = file.path("DEG_lists/adjusted_model/all_genes",
                     paste0(name, "_adjusted_allgenes.csv"))
  )
}
#save significant deg lists as .csv files
for(name in names(deg_lists)){
  sig <- subset(
    deg_lists[[name]],
    padj < 0.05 &
      abs(log2FoldChange) > 1
  )
  
  write.csv(
    as.data.frame(sig),
    file = file.path("DEG_lists/adjusted_model/sig_FC1",
                     paste0(name, "_adjusted_sigFC1.csv"))
  )
}

#save DESeq2 results object
saveRDS(dds, "deseq2_adjustedmodel_dds.rds")

#obtain results for co-variates
res_immune <- results(dds, name="immune.ESTIMATE")
res_stromal <- results(dds, name="stromal.ESTIMATE")
#export DEG list as a .csv file
write.csv(as.data.frame(res_stromal), file="stromal_results.csv")
write.csv(as.data.frame(res_immune), file="immune_results.csv")
#only significant DEGs (padj<0.05)
res_stromal.sig <- subset(res_stromal, padj < 0.05)
write.csv(as.data.frame(res_stromal.sig), file="stromal_deg_sig.csv")
res_immune.sig <- subset(res_immune, padj < 0.05)
write.csv(as.data.frame(res_immune.sig), file="immune_deg_sig.csv")