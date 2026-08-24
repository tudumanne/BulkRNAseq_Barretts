#Figure 2B
#This script contains the R code used to run decoupleR pathway analysis based on the following tutorial.
#https://saezlab.github.io/decoupleR/articles/pw_bk.html

#Input files: DEG lists with gene names and log normalized counts 

#load libraries
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(OmnipathR)

#load the DEG list 
deg_list = read.csv("deg_list_genename.csv")
rownames(deg_list) <- deg_list[,1]
deg_list <- deg_list[,-1]
deg_list <- na.omit(deg_list)
deg = deg_list %>% select(stat)
head(deg)

#load counts
counts = read.csv("normalized_counts_log.csv", header=TRUE)
#subset counts as required here based on sample ID
#N=NSq, B=NDBE, L=LGD, H=HGD, C=EAC
subset_counts = counts %>% select(gene_name | starts_with("N") |starts_with("B") |starts_with("L") |starts_with("H") | starts_with("C"))
rownames(subset_counts) <- subset_counts[,1]
subset_counts <- subset_counts[,-1]
head(subset_counts)

#pathway analysis
net <- decoupleR::get_progeny(organism = 'human', 
                              top = 500)

net

sample_acts <- decoupleR::run_mlm(mat = subset_counts, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'weight', 
                                  minsize = 5)
sample_acts

#transform to wide matrix
sample_acts_mat <- sample_acts %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

#scale per feature
sample_acts_mat <- scale(sample_acts_mat)

#color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05,2, length.out = floor(100 / 2)))

#plot the heatmap
pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 20,
                   cellheight = 20,
                   treeheight_row = 20,
                   treeheight_col = 20)

#run mlm
contrast_acts <- decoupleR::run_mlm(mat  =deg, 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor = 'weight', 
                                    minsize = 5)

contrast_acts

#plot the scores (Figure 2B in the manuscript)
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("Pathways")

p

#plot individual pathway of interest
pathway <- 'TNFa'

df <- net %>%
  dplyr::filter(source == pathway) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, 
                color = "3") %>%
  tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))

df <- df[inter, ]

df['stat'] <- deg[inter, ]

df <- df %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & stat > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & stat < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & stat > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & stat < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = df, 
                     mapping = ggplot2::aes(x = weight, 
                                            y = stat, 
                                            color = color)) + 
  ggplot2::geom_point(size = 2.5, 
                      color = "black") + 
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
  ggplot2::ggtitle(pathway)

p