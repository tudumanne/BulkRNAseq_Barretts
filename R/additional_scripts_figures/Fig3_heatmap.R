#Manuscript figure 3A
#R script for generating heatmap based on log2 normalized counts
#Reference: 
#Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.
#https://github.com/kevinblighe/E-MTAB-6141

require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)

library("dplyr")

deseq_counts = read.csv("normalized_counts_log.csv", header = TRUE)

gene_list = read.csv("genelist_heatmap_all.csv", header = TRUE)

meta_data = read.csv("metadata_heatmap.csv", header = TRUE)

deseq_counts = deseq_counts %>% 
  rename("gene_name"="V3","gene_id"="V2")

deseq_counts_selected = gene_list %>% inner_join(deseq_counts, by = "gene_id")

deseq_counts_selected2 = deseq_counts_selected[,!names(deseq_counts_selected) %in% c("gene_name.x","gene_name.y")]

rownames(deseq_counts_selected2) <- deseq_counts_selected2[,1]
deseq_counts_selected3 <- deseq_counts_selected2[,-1]

rownames(meta_data) <- meta_data[,1]

#dim(deseq_counts_selected3)
#dim(meta_data)
#all(rownames(meta_data) == colnames(deseq_counts_selected3))

meta_data$Classification = factor(meta_data$Classification, levels = c("NSq", "NDBE", "LGD", "HGD", "EAC"))

# Identify sample columns only
sample_cols <- rownames(meta_data)

expr_mat <- deseq_counts_selected3[, sample_cols]

# sanity check
stopifnot(colnames(expr_mat) == rownames(meta_data))

gene_ann_df <- deseq_counts_selected3 %>%
  dplyr::select(group, gene_type)

rownames(gene_ann_df) <- rownames(deseq_counts_selected3)

stopifnot(rownames(gene_ann_df) == rownames(expr_mat))

gene_ann_df$group = factor(gene_ann_df$group, levels = c("upregulated", "downregulated"))
gene_ann_df$gene_type = factor(gene_ann_df$gene_type, levels = c("protein_coding", "ncRNA"))
#max(deseq_counts_selected3)
#heat <- t(scale(t(deseq_counts_selected3)))

#myCol <- colorRampPalette(c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081'))(1000)
myCol <- colorRampPalette(c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081'))(1000)
myBreaks <- seq(0, 19, length.out = 1000)


colAnn <- HeatmapAnnotation(
  Classification = meta_data$Classification,
  col = list(
    Classification = c(
      "NSq"="#984ea3",
      "IM"="#377eb8",
      "LGD"="#4daf4a",
      "HGD"="#a65628",
      "Cancer"="#e41a1c"
    )
  )
)

rowAnn <- rowAnnotation(
  df = gene_ann_df,
  col = list(
    group = c(
      "upregulated" = "#EF8A62",
      "downregulated" = "#76A7CB"
    ),
    gene_type = c(
      "protein_coding" = "#984ea3",
      "ncRNA" = "#999999"
    )
  ))

hmap <- Heatmap(
  expr_mat,
  name = "log2 normalised counts",
  col = colorRamp2(myBreaks, myCol),
  top_annotation = colAnn,
  right_annotation = rowAnn,
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)

draw(hmap,
     heatmap_legend_side = 'bottom',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')
#dev.off()
