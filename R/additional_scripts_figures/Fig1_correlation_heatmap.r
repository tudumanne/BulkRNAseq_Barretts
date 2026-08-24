#Manuscript figure 1B
#R script for generating a correlation heatmap 

library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

dds <- readRDS("deseq2_simplemodel_dds.rds")
vsd <- vst(dds)

dds$Classification <- factor(
  dds$Classification,
  levels = c("NSq", "NDBE", "LGD", "HGD", "EAC")
)

top_var_genes <- head(
  order(rowVars(assay(vsd)), decreasing = TRUE),
  500
)
mat <- assay(vsd)[top_var_genes, ]

cor_mat <- cor(mat, method = "pearson")
cor_mat[cor_mat == 1] <- NA

sample_order <- order(vsd$Classification)
cor_mat <- cor_mat[sample_order, sample_order]

stage_cols <- c(
  "NSq"="#984ea3", "NDBE"="#377eb8", "LGD"="#4daf4a", "HGD"="#a65628", "EAC"="#e41a1c"
)

ha_top <- HeatmapAnnotation(
  Stage = vsd$Classification[sample_order],
  col = list(Stage = stage_cols),
  annotation_name_side = "left"
)

ha_left <- rowAnnotation(
  Stage = vsd$Classification[sample_order],
  col = list(Stage = stage_cols),
  annotation_name_side = "top"
)


col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#2166ac", "white", "#b2182b")
)

ht <- Heatmap(
  cor_mat,
  name = "Pearson\nr",
  col = col_fun,
  top_annotation = ha_top,
  left_annotation = ha_left,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = vsd$Classification[sample_order],
  column_split = vsd$Classification[sample_order],
  #clustering_distance_rows = function(x) as.dist(1 - x),
  #clustering_distance_columns = function(x) as.dist(1 - x),
  heatmap_legend_param = list(
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1")))

draw(ht)

# min -0.658596119
# max 0.97548174
