#Manuscript supplementary figure S1D-G
#R script for generating correlation scatter plots

library(ggplot2)
library(dplyr)
library(ggpubr)

deg = read.csv("NDBE_vs_NSq_deg.csv", header=TRUE)
deg <- deg[,!names(deg) %in% c("baseMean","lfcSE","stat","pvalue")]
deg <- deg %>% 
  rename("simple.log2FC"="log2FoldChange", "simple.padj"="padj")

cov = read.csv("adjusted_NDBE_vs_NSq_deg.csv", header=TRUE)
cov <- cov[,!names(cov) %in% c("baseMean","lfcSE","stat","pvalue")]
cov <- cov %>% 
  rename("adjusted.log2FC"="log2FoldChange", "adjusted.padj"="padj")

overlap = deg %>% inner_join(cov, by = "X")

overlap = na.omit(overlap)

overlap = overlap %>% mutate(p.val.cutoff =
                     case_when((adjusted.padj < 0.05 & simple.padj < 0.05) ~ "both",
                               simple.padj < 0.05 ~ "simple-only", 
                               adjusted.padj < 0.05 ~ "adjusted-only",
                                    TRUE ~ "neither"))

ggplot(overlap, aes(x=simple.log2FC, y=adjusted.log2FC, add="reg.line"))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(color=p.val.cutoff, fill=p.val.cutoff), size=2, alpha=0.5) +
  scale_color_manual(values=c("both"="brown", "simple-only"="blue", "adjusted-only"="purple", "neither"="grey"))+
  geom_smooth(method='lm', se = TRUE, linewidth = 1)+
  xlim(c(-10,10))+
  ylim(c(-10,10))+
  xlab("RNAseq-simple - log2FC")+
  ylab("RNAseq-adjusted - log2FC")+
  stat_cor(method = "pearson", cor.coef.name = "r", vjust = 1, size = 4) +
  #stat_regline_equation(label.y = 4.2)+
  theme_bw()+
  theme(text = element_text(colour = "black"), plot.title=element_text(face="bold"), plot.subtitle=element_text(face="bold"), plot.caption =element_text(face="bold"),
        axis.title.x=element_text(face="bold", size=12, colour = "black"), axis.text.x = element_text(size=12, colour = "black"),
        axis.title.y=element_text(face="bold", size=12, colour = "black"), axis.text.y = element_text(size=12, colour = "black"), 
        legend.text = element_text(size=12, colour = "black"), legend.title = element_text(face = "bold", size=12, colour = "black"))

#double check if the new column added is correct