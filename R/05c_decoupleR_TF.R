#Figure 2C
#This script contains the R code used to run decoupleR transcription factor (TF) analysis based on the following tutorial.
#https://saezlab.github.io/decoupleR/articles/tf_bk.html

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
deg = read.csv("deg_list_genename.csv")
rownames(deg) <- deg[,1]
deg <- deg[,-1]
deg <- na.omit(deg)
head(deg)

counts = read.csv("normalized_counts_log.csv", header=TRUE)
#subset counts as required here based on sample ID
#N=NSq, B=NDBE, L=LGD, H=HGD, C=EAC
subset_counts = counts %>% select(gene_name | starts_with("N") | starts_with("B") |starts_with("L") |starts_with("H") | starts_with("C"))
rownames(subset_counts) <- subset_counts[,1]
subset_counts <- subset_counts[,-1]
head(subset_counts)

#TF analysis
net <- get_collectri(organism='human', split_complexes=FALSE)
net

sample_acts <- run_ulm(mat=subset_counts, net=net, .source='source', .target='target',
                       .mor='mor', minsize = 5)
sample_acts

n_tfs <- 25

#transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

#get top tfs with more variable means across clusters
tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
sample_acts_mat <- sample_acts_mat[,tfs]

#scale per sample
sample_acts_mat <- scale(sample_acts_mat)

#choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

#plot the heatmap
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)

#run ulm
contrast_acts <- run_ulm(mat=deg[, 'stat', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)
contrast_acts

#filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

#plot the scores (Figure 2C in the manuscript)
ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("TFs")

#plot individual TF
tf <- 'HNF1A'

df <- net %>%
  filter(source == tf) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('logfc', 'stat', 'p_value')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & stat > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & stat < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & stat < 0, '1', color))

ggplot(df, aes(x = logfc, y = -log10(p_value), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(tf)
