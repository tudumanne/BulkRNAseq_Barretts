#Manuscript supplementary figure S4
#R script for generating violin plots showing log2 normalized counts

counts = read.csv("normalized_counts_log.csv", header=TRUE)

counts2 = counts[,-1]

library(data.table) #load the libraries in this order, same function in another package
counts3 = transpose(counts2, keep.names="sample")

library(janitor)
counts3 = counts3 %>%
  row_to_names(row_number=1)

colnames(counts3)[1] ="sample"

#counts3$sample<-gsub("X","",as.character(counts3$sample))


genelist = read.csv("upregulated_shortlist.csv", header=TRUE)


library(tidyverse)
genelist_counts = counts3 %>% select(c(sample, genelist$gene_name))

sample_categories = read.csv("metadata_heatmap.csv", header=TRUE)

sample_categories = sample_categories %>% 
  rename("sample"="X")

genelist_counts2 = genelist_counts %>% inner_join(sample_categories, by = "sample")

library(ggplot2)

genelist_counts2$Classification = factor(genelist_counts2$Classification, levels=c("NSq", "NDBE", "LGD", "HGD", "EAC"))

genelist_counts2[, c(2:14)] <- sapply(genelist_counts2[, c(2:14)], as.numeric)

loop.cols <- names(genelist_counts2)[2:14]
for (col in loop.cols) {
  p <- ggplot(genelist_counts2, aes(x=Classification, y=!!sym(col), fill=Classification)) + 
    geom_violin(alpha=0.3, linewidth=0.5, trim=FALSE, aes(fill=Classification))+
    scale_fill_manual(values=c("NSq"="#984ea3", "NDBE"="#377eb8", "LGD"="#4daf4a", "HGD"="#a65628", "EAC"="#e41a1c"))+
    labs(title=col, x="Classification", y = "log2 normalised counts")+
    geom_boxplot(width=0.1)+
    theme_bw()
    print(p)
    ggsave(filename=paste0(col, ".pdf"), p, width = 6, height = 4)
}
dev.off()


