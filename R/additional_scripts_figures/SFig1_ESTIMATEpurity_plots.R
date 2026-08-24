#Manuscript supplementary figure S1A-C
#R script for generating box plots

library(ggplot2)

data = read.csv("ESTIMATE_purity_scores_85.csv", header = TRUE)

data$classification = factor(data$classification, levels=c("NSq", "NDBE", "LGD", "HGD", "EAC"))

ggplot(data, aes(classification, immune.ESTIMATE, fill=classification))+ 
  geom_boxplot(alpha=0.3, linewidth=0.5)+
  scale_fill_manual(values=c("NSq"="#984ea3", "NDBE"="#377eb8", "LGD"="#4daf4a", "HGD"="#a65628", "EAC"="#e41a1c"))+
  theme_bw()+
  scale_y_continuous(limits = c(-900, 1250), breaks = seq(-900, 1250, 500))

ggplot(data, aes(classification, stromal.ESTIMATE, fill=classification))+ 
  geom_boxplot(alpha=0.3, linewidth=0.5)+
  scale_fill_manual(values=c("NSq"="#984ea3", "NDBE"="#377eb8", "LGD"="#4daf4a", "HGD"="#a65628", "EAC"="#e41a1c"))+
  theme_bw()+
  scale_y_continuous(limits = c(-2250, 1150), breaks = seq(-2250, 1150, 1000))


ggplot(data, aes(classification, purity, fill=classification))+ 
  geom_boxplot(alpha=0.3, linewidth=0.5)+
  scale_fill_manual(values=c("NSq"="#984ea3", "NDBE"="#377eb8", "LGD"="#4daf4a", "HGD"="#a65628", "EAC"="#e41a1c"))+
  theme_bw()+
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2))





