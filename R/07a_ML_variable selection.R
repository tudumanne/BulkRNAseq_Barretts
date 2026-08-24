#This script contains the R code used for variable selection (Machine learning based ranking of biomarkers)
#Input files: sample metadata, log normalized count data from DESeq2 and TPM counts from Salmon

##################################
# Authro: Dana Pascovici

# Variable selection using
#  xboost 
#  lasso 
# check accuracy and output results to csv
##################################
rm(list=ls())

library(openxlsx)
library(MASS)
library(caret)
library(pROC)
# various "classic" ML packages
library(xgboost)
library(glmnet)  # for lasso regression
library(nnet)

# For setting up, input:
# data matrix (variables as rows, samples as columns)

tpm = read.csv("data/counts_tpm.csv")
rownames(tpm) = tpm[,1]
tpm = tpm[,-1]
tpm.gp =  substring(colnames(tpm), 1,1)
tpm.ag = aggregate(t(tpm), by=list(Group = tpm.gp), FUN=median)

rownames(tpm.ag) = tpm.ag[,1]
tpm.ag = tpm.ag[,-1]
max.tpm = apply(tpm.ag, 2, FUN=max)
filter.idx = max.tpm > 5
filter.vals = colnames(tpm.ag)[filter.idx]
filter.vals = gsub("\\..*", "", filter.vals)

dat = read.csv(file.path("data", "discovery_biomarkers_log_counts.csv"))

keep.idx = dat[,1] %in% filter.vals
dat = dat[keep.idx,]

dd = dat[,-2]        # data with variables names in column 1
markers = dat[,1:2]  # any info on the variables
pd = data.frame(Sample=colnames(dat)[-c(1:2)])  # any phenotypic data
pd$Condition = substring(pd$Sample, 1,1)
table(pd$Condition)

# change to the appropriate comparison here; N+B+L vs H+C OR  N+B vs L+H+C
# N=NSq, B=NDBE, L=LGD, H=HGD and C=EAC
pd$Status = 1
pd$Status[pd$Condition %in% c("N", "B")] = 0

mydat = data.frame(as.numeric(as.factor(pd$Status))-1, t(dd[,-1]))
colnames(mydat) = c("Status", paste0(dd[,1]))
mydat.all = mydat

#############################
# input to ML algorithms: 
# data (samples x variables)
# first column: status
#############################

# variable selection with lasso / elastic net
# Description link:
# http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/
# x =	input matrix, of dimension nobs x nvars; each row is an observation vector.
# y = response variable. 
x = as.matrix(mydat[,-1])
y = mydat$Status
y_xg = as.factor(mydat$Status)

#set seed
set.seed(12345)

lasso.list = list()

for (iii in 1:100){
  
  rand = sample(6, nrow(x), replace=T)
  cv.lasso <- cv.glmnet(x[rand!=1,], y[rand!= 1], alpha = 1, family = "binomial") # see names for various values, including lambda values
  
  # plot(cv.lasso)
  ccc = coef(cv.lasso, cv.lasso$lambda.min)
  ccc = coef(cv.lasso, cv.lasso$lambda.1se)
  lasso.list[[iii]] = rownames(ccc)[which(ccc != 0)]
}

# variable selection with xgboost
xgboost.list = list()

# variable selection with xboost 
# Vignette link:
# https://cran.r-project.org/web/packages/xgboost/vignettes/discoverYourData.html
for (iii in 1:100){
  
rand = sample(6, nrow(x), replace=T)
bst <- xgboost(data =x[rand!=1,], label = y_xg[rand!= 1], max_depth = 4,
               eta = 1, nthread = 2, nrounds = 10,objective = "binary:logistic")

importance <- xgb.importance(feature_names = colnames(x[rand!=1,]), model = bst)
xgboost.list[[iii]] = importance$Feature

}


#select top 20 
selectN = 20

##############################################
# Output: candidates, importance, model, etc
##############################################

lasso.table = sort(table(unlist(lasso.list)), decreasing=TRUE)[-1][1:selectN]   # drop intercept!
xgboost.table = sort(table(unlist(xgboost.list)), decreasing=TRUE)[1:selectN]

importance.summaries = data.frame(lasso.table, xgboost.table)
colnames(importance.summaries)[c(1,3)] = c("Lasso", "XGboost")
importance.summaries

write.csv(importance.summaries, file=file.path("results", "VarSelectionResults_LHCvsNB.csv"))
