# Load library
library("dplyr")

counts = read.csv("normalized_counts_untransformed_dataset1.csv", header=TRUE)

counts <- counts %>% 
  rename("gene_id"="V2", "gene_name"="V3")

subset_counts_nsq = counts %>% select(gene_id |gene_name | starts_with("N"))
subset_counts_nsq$nsq_avg_counts = rowMeans(subset_counts_nsq[c(3:12)])
subset_counts_nsq$nsq_presence = rowSums(subset_counts_nsq[c(3:12)] > 0)
counts_nsq = subset_counts_nsq %>%
  select(gene_id, gene_name, nsq_avg_counts, nsq_presence)

subset_counts_im = counts %>% select(gene_id | starts_with("B"))
subset_counts_im$im_avg_counts = rowMeans(subset_counts_im[c(2:30)])
subset_counts_im$im_presence = rowSums(subset_counts_im[c(2:30)] > 0)
counts_im = subset_counts_im %>%
  select(gene_id, im_avg_counts, im_presence)

subset_counts_lgd = counts %>% select(gene_id | starts_with("L"))
subset_counts_lgd$lgd_avg_counts = rowMeans(subset_counts_lgd[c(2:22)])
subset_counts_lgd$lgd_presence = rowSums(subset_counts_lgd[c(2:22)] > 0)
counts_lgd = subset_counts_lgd %>%
  select(gene_id, lgd_avg_counts, lgd_presence)

subset_counts_hgd = counts %>% select(gene_id | starts_with("H"))
subset_counts_hgd$hgd_avg_counts = rowMeans(subset_counts_hgd[c(2:12)])
subset_counts_hgd$hgd_presence = rowSums(subset_counts_hgd[c(2:12)] > 0)
counts_hgd = subset_counts_hgd %>%
  select(gene_id, hgd_avg_counts, hgd_presence)

subset_counts_cancer = counts %>% select(gene_id | starts_with("C"))
subset_counts_cancer$cancer_avg_counts = rowMeans(subset_counts_cancer[c(2:15)])
subset_counts_cancer$cancer_presence = rowSums(subset_counts_cancer[c(2:15)] > 0)
counts_cancer = subset_counts_cancer %>%
  select(gene_id, cancer_avg_counts, cancer_presence)

subset_counts_imdyscancer = counts %>% select(gene_id | starts_with("B") | starts_with("L")| starts_with("H")| starts_with("C"))
subset_counts_imdyscancer$imdyscancer_avg_counts = rowMeans(subset_counts_imdyscancer[c(2:76)])
subset_counts_imdyscancer$imdyscancer_presence = rowSums(subset_counts_imdyscancer[c(2:76)] > 0)
counts_imdyscancer = subset_counts_imdyscancer %>%
  select(gene_id, imdyscancer_avg_counts, imdyscancer_presence)

subset_counts_dyscancer = counts %>% select(gene_id | starts_with("L")| starts_with("H")| starts_with("C"))
subset_counts_dyscancer$dyscancer_avg_counts = rowMeans(subset_counts_dyscancer[c(2:47)])
subset_counts_dyscancer$dyscancer_presence = rowSums(subset_counts_dyscancer[c(2:47)] > 0)
counts_dyscancer = subset_counts_dyscancer %>%
  select(gene_id, dyscancer_avg_counts, dyscancer_presence)

subset_counts_hgdcancer = counts %>% select(gene_id | starts_with("H")| starts_with("C"))
subset_counts_hgdcancer$hgdcancer_avg_counts = rowMeans(subset_counts_hgdcancer[c(2:26)])
subset_counts_hgdcancer$hgdcancer_presence = rowSums(subset_counts_hgdcancer[c(2:26)] > 0)
counts_hgdcancer = subset_counts_hgdcancer %>%
  select(gene_id, hgdcancer_avg_counts, hgdcancer_presence)

subset_counts_imlgd = counts %>% select(gene_id | starts_with("B") | starts_with("L"))
subset_counts_imlgd$imlgd_avg_counts = rowMeans(subset_counts_imlgd[c(2:51)])
subset_counts_imlgd$imlgd_presence = rowSums(subset_counts_imlgd[c(2:51)] > 0)
counts_imlgd = subset_counts_imlgd %>%
  select(gene_id, imlgd_avg_counts, imlgd_presence)

subset_counts_dys = counts %>% select(gene_id | starts_with("L")| starts_with("H"))
subset_counts_dys$dys_avg_counts = rowMeans(subset_counts_dys[c(2:33)])
subset_counts_dys$dys_presence = rowSums(subset_counts_dys[c(2:33)] > 0)
counts_dys = subset_counts_dys %>%
  select(gene_id, dys_avg_counts, dys_presence)

counts1 = counts_nsq %>% inner_join(counts_im, by = "gene_id")
counts2 = counts1 %>% inner_join(counts_lgd, by = "gene_id")
counts3 = counts2 %>% inner_join(counts_hgd, by = "gene_id")
counts4 = counts3 %>% inner_join(counts_cancer, by = "gene_id")
counts5 = counts4 %>% inner_join(counts_imdyscancer, by = "gene_id")
counts6 = counts5 %>% inner_join(counts_dyscancer, by = "gene_id")
counts7 = counts6 %>% inner_join(counts_hgdcancer, by = "gene_id")
counts8 = counts7 %>% inner_join(counts_imlgd, by = "gene_id")
counts_summary = counts8 %>% inner_join(counts_dys, by = "gene_id")


write.csv(counts_summary, "counts_summary.csv")
save(counts_summary, file="counts_summary.RData")

salmon_tpm = read.csv("salmon_counts.csv", header=TRUE)

salmon_tpm$nsq_median_tpm = apply(salmon_tpm[,grep("N",colnames(salmon_tpm))], 1, median, na.rm=T)
salmon_tpm$im_median_tpm = apply(salmon_tpm[,grep("B",colnames(salmon_tpm))], 1, median, na.rm=T)
salmon_tpm$lgd_median_tpm = apply(salmon_tpm[,grep("L",colnames(salmon_tpm))], 1, median, na.rm=T)
salmon_tpm$hgd_median_tpm = apply(salmon_tpm[,grep("H",colnames(salmon_tpm))], 1, median, na.rm=T)
salmon_tpm$cancer_median_tpm = apply(salmon_tpm[,grep("C",colnames(salmon_tpm))], 1, median, na.rm=T)
imdyscancer = c("B","L", "H", "C")
salmon_tpm$imdyscancer_median_tpm = apply(salmon_tpm[,grep(paste(imdyscancer,collapse="|"),colnames(salmon_tpm))], 1, median, na.rm=T)
dyscancer = c("L","H","C")
salmon_tpm$dyscancer_median_tpm = apply(salmon_tpm[,grep(paste(dyscancer,collapse="|"),colnames(salmon_tpm))], 1, median, na.rm=T)
hgdcancer = c("H","C")
salmon_tpm$hgdcancer_median_tpm = apply(salmon_tpm[,grep(paste(hgdcancer,collapse="|"),colnames(salmon_tpm))], 1, median, na.rm=T)
lgdim = c("B","L")
salmon_tpm$imlgd_median_tpm = apply(salmon_tpm[,grep(paste(lgdim,collapse="|"),colnames(salmon_tpm))], 1, median, na.rm=T)
dys = c("L","H")
salmon_tpm$dys_median_tpm = apply(salmon_tpm[,grep(paste(dys,collapse="|"),colnames(salmon_tpm))], 1, median, na.rm=T)

tpm_summary = salmon_tpm %>%
  select(gene_id, nsq_median_tpm, im_median_tpm, lgd_median_tpm, hgd_median_tpm, cancer_median_tpm,
         imlgd_median_tpm, imdyscancer_median_tpm, dyscancer_median_tpm, hgdcancer_median_tpm, dys_median_tpm)

write.csv(tpm_summary, "tpm_summary.csv")
save(tpm_summary, file="tpm_summary.RData")


deg = read.csv("im_nsq_deg_padj0.05_fc1.csv", header = TRUE)
deg_cov = read.csv("cov_im_nsq_deg_padj0.05_fc1.csv", header = TRUE)

deg_cov <- deg_cov %>% 
  rename("log2FoldChange.cov"="log2FoldChange","padj.cov"="padj")

deg2 <- deg[,!names(deg) %in% c("baseMean","lfcSE","stat","pvalue")]
deg_cov2 <- deg_cov[,!names(deg_cov) %in% c("baseMean","lfcSE","stat","pvalue")]


deg_overlap = deg2 %>% inner_join(deg_cov2, by = "X")
deg_overlap <- deg_overlap %>% 
  rename("gene_id"="X")

genes = read.csv("tx2gene.gencode.v36.pa_edited.csv", header=TRUE)

deg_overlap2 = deg_overlap %>% inner_join(genes, by = "gene_id")

load("counts_summary.RData")
load("tpm_summary.RData")

deg_overlap2 <- deg_overlap2 %>% 
  left_join(select(counts_summary, gene_name, nsq_avg_counts, nsq_presence, im_avg_counts, im_presence), by = "gene_name")

deg_overlap2$imnsq_ratio_counts = deg_overlap2$im_avg_counts/deg_overlap2$nsq_avg_counts

deg_overlap2 <- deg_overlap2 %>% 
  left_join(select(tpm_summary, gene_id, nsq_median_tpm, im_median_tpm), by = "gene_id")

imnsq_deg_overlap = select(deg_overlap2, gene_id, gene_name, log2FoldChange, padj, 
                                   log2FoldChange.cov, padj.cov, nsq_avg_counts, im_avg_counts, imnsq_ratio_counts,
                                   nsq_median_tpm, im_median_tpm, nsq_presence, im_presence)

imnsq_subset_up_shortlist = imnsq_deg_overlap %>% filter(log2FoldChange >= 2 & imnsq_ratio_counts >= 5 & im_median_tpm >= 5)
imnsq_subset_down_shortlist = imnsq_deg_overlap %>% filter(log2FoldChange <= -2 & imnsq_ratio_counts <= 0.2 & im_median_tpm == 0)

write.csv(imnsq_subset_up_shortlist, "imnsq_upregulated_shortlist.csv")
write.csv(imnsq_subset_down_shortlist, "imnsq_downregulated_shortlist.csv")

save(imnsq_deg_overlap, imnsq_subset_up_shortlist, imnsq_subset_down_shortlist, file="im.vs.nsq.RData")



