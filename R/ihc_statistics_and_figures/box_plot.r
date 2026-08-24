library(dplyr)
library(tidyverse)
library(ggplot2)

h_score = read.csv("all_scores_edited.csv", header=TRUE)
metadata = read.csv("all_metadata.csv", header=TRUE)
h_score_merged = h_score %>% inner_join(metadata, by = "Tissue_ID")

h_score_merged$Type = factor(h_score_merged$Type, levels=c("Normal","Malignant"))
h_score_merged$TP53 = as.numeric(h_score_merged$TP53)

mw <- wilcox.test(TP53 ~ Type, data = h_score_merged)
mw
p_val <- mw$p.value
W_val <- mw$statistic

n1 <- sum(h_score_merged$Type == levels(h_score_merged$Type)[1])
n2 <- sum(h_score_merged$Type == levels(h_score_merged$Type)[2])

r_rb <- (2 * W_val) / (n1 * n2) - 1
r_rb <- round(r_rb, 2)
r_rb

stars <- ifelse(
  p_val < 0.0001, "****",
  ifelse(
    p_val < 0.001, "***",
    ifelse(
      p_val < 0.01, "**",
      ifelse(p_val < 0.05, "*", "ns")
    )
  )
)

annot_h_score_merged <- data.frame(
  x_start = 1,
  x_end   = 2,
  y       = 320,
  stars   = stars
)


#portrait 4x5
ggplot(h_score_merged, aes(x=Type, y=TP53, fill=Type)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  #geom_violin(alpha=0.3, linewidth=0.5, trim=FALSE, aes(fill=Type))+
  scale_fill_manual(values=c("Normal"="#66c2a4", "Malignant"="#fcae91"))+
  labs(title="TP53",x="Tissue category", y = "H-score")+
  #geom_boxplot(width=0.1)+
  theme_bw()+
  ylim(-20,350)+
#geom_jitter(shape=20, position=position_jitter(0.2))+
  
  # horizontal bracket
  geom_segment(
    data = annot_h_score_merged,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  
  # vertical ticks
  geom_segment(
    data = annot_h_score_merged,
    aes(x = x_start, xend = x_start, y = y, yend = y * 0.97),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = annot_h_score_merged,
    aes(x = x_end, xend = x_end, y = y, yend = y * 0.97),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  
  # stars
  geom_text(
    data = annot_h_score_merged,
    aes(x = 1.5, y = y * 1.03, label = stars),
    inherit.aes = FALSE,
    size = 4
  )

