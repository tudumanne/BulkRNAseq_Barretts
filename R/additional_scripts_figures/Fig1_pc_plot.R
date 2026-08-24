#Manuscript figure 1A
#R script for generating PCA plot

library(DESeq2)
library(ggplot2)

#transformation of data for visualization
vsd <- vst(dds)

#define classification
dds$Classification <- factor(
  dds$Classification,
  levels = c("NSq", "NDBE", "LGD", "HGD", "EAC")
)

#PCA
pcaData <- plotPCA(vsd, intgroup=c("Classification"), returnData=TRUE, ntop=500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(data = pcaData, aes(PC1, PC2)) +
  geom_point(aes(color=Classification, fill=Classification), size=3, alpha=0.8) +
  scale_color_manual(values=c("NSq"="#984ea3", "NDBE"="#377eb8", "LGD"="#4daf4a", "HGD"="#a65628", "EAC"="#e41a1c"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  #geom_label_repel(aes(label = rownames(pcaData), size = 1)) +
  theme_bw()+
  theme(text = element_text(colour = "black"), plot.title=element_text(face="bold"), plot.subtitle=element_text(face="bold"), plot.caption =element_text(face="bold"),
        axis.title.x=element_text(face="bold", size=12, colour = "black"), axis.text.x = element_text(face = "bold", size=12, colour = "black"),
        axis.title.y=element_text(face="bold", size=12, colour = "black"), axis.text.y = element_text(face = "bold", size=12, colour = "black"), 
        legend.text = element_text(face = "bold", size=12, colour = "black"), legend.title = element_text(face = "bold", size=12, colour = "black") )
p