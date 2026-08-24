#This script contains the R code used for tissue purity estimation using ESTIMATE algorithm, 
#implemented via tidyestimate R package.
#Input files: normalized counts from DESeq2 analysis (simple model)

#References:
#https://github.com/KaiAragaki/tidyestimate

library(tidyestimate)
library(dplyr)

data = read.table("normalized_counts.txt", header = TRUE, sep="\t")

dim(data)
head(data)[,1:5]

data2 = data %>% distinct(V3, .keep_all=TRUE)
data3 = data.frame(data2, row.names = 2)
data4 = data3[ ,-c(1)]

scores = data4 |>
    filter_common_genes(id="hgnc_symbol", tell_missing=FALSE, find_alias=TRUE) |>
    estimate_score(is_affymetrix=TRUE)

#is_affymetrix=FALSE will not generate purity scores
#https://rdrr.io/cran/tidyestimate/man/estimate_score.html

scores |>
    plot_purity(is_affymetrix=TRUE)

write.csv(scores, "estimate_scores.csv")
