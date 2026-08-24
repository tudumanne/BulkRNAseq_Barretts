#Manuscript supplementary figure 2A-B
#R script for generating UpSet plots

library(ComplexUpset)
library(dplyr)
library(tidyr)
library(ggplot2)

ndbe_nsq = read.csv("NDBE_Nsq_up.csv", header=TRUE)
lgd_nsq = read.csv("LGD_Nsq_up.csv", header=TRUE)
hgd_nsq = read.csv("HGD_Nsq_up.csv", header=TRUE)
eac_nsq = read.csv("EAC_Nsq_up.csv", header=TRUE)

set_a = ndbe_nsq$X
set_b = lgd_nsq$X
set_c = hgd_nsq$X
set_d = eac_nsq$X

all_items <- unique(c(set_a, set_b, set_c, set_d))

df <- tibble(
  item = all_items,
  a = item %in% set_a,
  b = item %in% set_b,
  c = item %in% set_c,
  d = item %in% set_d
)

p <- upset(
  df,
  intersect = c("a", "b", "c", "d"),
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  set_sizes = upset_set_size() + 
    geom_text(
      aes(label = after_stat(count)), 
      stat = 'count',
      hjust = 1.1,          
      color = "black",      
      size = 3.5            
    ) +
    theme(
      plot.margin = margin(t = 0, r = 0, b = 0, l = 2, unit = "cm")
    ),
  width_ratio = 0.2
) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p + theme(plot.margin = unit(c(0, 0, 0, 1), "cm"))
