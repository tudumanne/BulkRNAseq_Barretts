#This script contains the R code used for ranking, feature selection and panel optimization (Machine learning based ranking of biomarkers)
#Input files: sample metadata, log normalized count data from both discovery and validation cohorts

##################################
# Author: Dana Pascovici
# 
# Use variable selection 
# generate ranking
# optimize a few combination panels
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
#################
# utility functions
###################

# CV functions from Venables and Ripley
CVtest = function(vec, fitfn, predfn, ...) {
  res = vec
  for (i in sort(unique(rand))) {
    # Uncomment line below if you want to see the CV folds as they are run
    # cat("Fold", i, "\n", sep="")
    learn = fitfn(subset = (rand != i), ...)
    res[rand==i] = predfn(learn)[rand==i]
  }
  res
}


con = function(...) {
  tab = table(...)
  print(tab)
  diag(tab) = 0
  cat("error rate = ", round(100*sum(tab)/length(list(...)[[1]]), 2), "%\n")
}


# input data set - discovery cohort
dat = read.csv(file.path("data", "discovery_biomarkers_log_counts.csv"))
rownames(dat) = dat[,1]

# input data set - test data: MAAG 2017 paper
dat.test = read.csv(file.path("data", "maag_validation_biomarkers_log_counts.csv"))
rownames(dat.test) = dat.test[,1]
test.metadata = read.csv(file.path("data", "maag_validation_metadata.csv"))
test.metadata$ConditionOrdered = factor(test.metadata$Classification, levels=c("NSq", "NDBE", "LGD", "EAC"))

########################
# training data processing
########################
markers = dat[,1:2]
rownames(markers) = markers[,1]

pd = data.frame(Sample=colnames(dat)[-c(1:2)])  # any phenotypic data
pd$Condition = substring(pd$Sample, 1,1)
pd$ConditionOrdered = factor(pd$Condition, levels=c("N", "B", "L", "H", "C"))
table(pd$Condition)

# change to the appropriate comparison here; N+B+L vs H+C OR  N+B vs L+H+C
# N=NSq, B=NDBE, L=LGD, H=HGD and C=EAC
pd$Status = 1
pd$Status[pd$Condition %in% c("N", "B")] = 0

# Data to use for classification
# Format: data frame with Status included, variables as columns
mydat = data.frame(as.numeric(as.factor(pd$Status))-1, t(dat[,-c(1:2)]))
colnames(mydat) = c("Status", dat$gene_id)
mydat.all = mydat

########################
# test data processing
# Maag dataset
########################

# change to the appropriate comparison here; N+B+L vs H+C OR  N+B vs L+H+C
# N=NSq, B=NDBE, L=LGD, and C=EAC
test.metadata$Status = 1
test.metadata$Status[test.metadata$Classification %in% c("NSq", "NDBE")] = 0

mydat.test = data.frame(Status = as.numeric(as.factor(test.metadata$Status))-1, t(dat.test[,-c(1:2)]))
colnames(mydat.test) = c("Status", dat.test[,1])

table(test.metadata$Classification)

# Use previously generated ranking
# Ranking generated using script "07a_ML_variable selection.R"

biomarkers = read.csv(file.path("results", "VarSelectionResults_LHCvsNB.csv"))

combine.ranks = function(biomarkers, marker.names)  {
  all.candidates=union(biomarkers[,2], biomarkers[,4])
  ranks = data.frame(lapply(c(2,4), FUN=function(i){match(all.candidates, biomarkers[,i], nomatch=21)}))
  rownames(ranks) = gsub("\\.", " ", all.candidates)
  candidates.ranked = data.frame(Candidate= all.candidates[order(apply(ranks, 1, FUN=mean))], MeanRank=floor(apply(ranks, 1, FUN=mean)[order(apply(ranks, 1, FUN=mean))]))
  biomarkers.ranked = candidates.ranked
  ranks.printed = 22-ranks
  rank.sum = apply(ranks.printed, 1, FUN=sum)
  colnames(ranks.printed) = c("LASSO", "XGboost")
  ranks.printed = data.frame(marker.names[rownames(ranks.printed),], ranks.printed, rank.sum)
  ranks.printed[order(rank.sum, decreasing=TRUE),]
}

biomarker.ranks = combine.ranks(biomarkers, markers)

AUC.train = sapply(biomarker.ranks[,1], FUN=function(x){ roc(pd$Status,  mydat[,x])$auc })

AUC.test = sapply(biomarker.ranks[,1], FUN=function(x){ roc(mydat.test$Status,  mydat.test[,x])$auc })

boxplot(data.frame(AUC.train, AUC.test), main="Top biomarkers - AUC train and test (Maag)", ylim=c(0.3,1))
abline(h=0.9, lwd=0.5, col="gray")

df = data.frame(biomarker.ranks, AUC.train, AUC.test)

write.csv(df, file=file.path("results", "TopRankedCandidates.csv"))

############################
# Panel optimization
############################

# previous ranking in 'biomarkers'
mydat = mydat.all

# consolidate LASSO and xgboost candidates
all.candidates=union(biomarkers[,2], biomarkers[,4])
ranks = data.frame(lapply(c(2,4), FUN=function(i){match(all.candidates, biomarkers[,i], nomatch=21)}))
rownames(ranks) = gsub("\\.", " ", all.candidates)
candidates.ranked = data.frame(Candidate= all.candidates[order(apply(ranks, 1, FUN=mean))], MeanRank=floor(apply(ranks, 1, FUN=mean)[order(apply(ranks, 1, FUN=mean))]))
AUC = sapply(candidates.ranked[,1], FUN=function(x){ roc(pd$Status,  mydat[,x])$auc })
biomarkers.ranked = candidates.ranked
ranks.printed = 22-ranks
rank.sum = apply(ranks.printed, 1, FUN=sum)

rownames(ranks.printed) = markers[rownames(ranks.printed),2]

# Ranking plot
#changed below 15 to 8
par(mar=c(3,8,2,2))
barplot(t(ranks.printed[order(rank.sum),]), horiz=TRUE, las=2, main="Variable ranking")
legend("bottomright", fill=grey.colors(2), legend=c("LASSO", "XGBoost"))

# top 9 ranked
top = markers[candidates.ranked[1:9,1],]

layout(matrix(1:9, nrow=3))
for (ii in 1:9) boxplot(mydat[,top[ii,1]] ~ pd$ConditionOrdered, main=top[ii,2]);

ranked.results = data.frame(candidates.ranked, markers[candidates.ranked[,1],], AUC)

#################################################
# Panel optimisation
# Approach 1: start with best performer
# try adding some more while performance improves
#################################################
# Reference:
# Thakur, Abhishek. Approaching (Almost) Any Machine Learning Problem (p. 160).
# Abhishek Thakur. 

#can start with top AUC or top ranked
#GoodFeatures = ranked.results[which.max(ranked.results$AUC), "Candidate"]
#Top ranked
GoodFeatures = ranked.results[1, "Candidate"]

# use starting point and whatever is available
Features = setdiff(all.candidates, GoodFeatures)

Group = as.numeric(as.factor(pd$Status)) - 1
# rand =  sample(6, length(Group), replace=T)  #1:length(Group)  #
rand = 1:length(Group)

mydat = mydat.all[,colnames(mydat.all) %in% GoodFeatures,drop=FALSE]
 
prob.logistic = CVtest( Group,
 function(x,...) glm(Group ~ ., mydat, , family=binomial(link="logit"), ...),
 function(obj, x) predict(obj, mydat, type='response') )


res.logistic = ifelse(prob.logistic>0.5, 1,0)
# optimise threshold
# res.logistic = as.factor(ifelse(prob.logistic>OpThresh(prob.logistic, as.factor(Group), Grid=seq(0.25,0.75,by=0.02), Criteria="acc", positive="1"), 1, 0))
con(res.logistic, Group)
Acc = sum(res.logistic != Group)/length(Group)

bestAUC = auc(Group, prob.logistic)

# Find which other features improves classification the most
# if added to the existing features up to this point
#initialize good features list  
#and best scores to keep track of both  
good_features = GoodFeatures
best_scores = rep(0, length(Features))
# best_scores[Features %in% GoodFeatures] = bestAUC

NN = 1 
while ( NN < 10 ) {

added = FALSE

cat(NN, "\n")
for (idx in 1:length(Features)) {
  
  ff = Features[idx]
  if (!( ff %in% good_features)) {
  current_features = c(good_features, ff)
  mydat = mydat.all[,colnames(mydat.all) %in% current_features,drop=FALSE]
  
  prob.logistic = CVtest( Group,
                          function(x,...) glm(Group ~ ., mydat, , family=binomial(link="logit"), ...),
                          function(obj, x) predict(obj, mydat, type='response') )
  
  
  res.logistic = ifelse(prob.logistic>0.5, 1,0)
  # optimise threshold
  # res.logistic = as.factor(ifelse(prob.logistic>OpThresh(prob.logistic, as.factor(Group), Grid=seq(0.25,0.75,by=0.02), Criteria="acc", positive="1"), 1, 0))
  # con(res.logistic, Group)
  
  f.auc = auc(Group, prob.logistic)
 
  best_scores[idx] = f.auc
  }
}
  
  
  if ( max(best_scores) >  bestAUC + 0.005 ) {  
    good_features = c(good_features, Features[which.max(best_scores)])
    Features = setdiff(Features,Features[which.max(best_scores)] )
    bestAUC = max(best_scores)
    added = TRUE
    NN = NN+1
  }
  

if (!added) { NN = 11 }

}

opt10 = c("gene1", "gene2", "gene3", "gene4", "gene5")

# panel performance in training data
Group = as.numeric(as.factor(pd$Status)) - 1
rand = 1:length(Group)
mydat = mydat.all[,colnames(mydat.all) %in% opt10,drop=FALSE]
# prob.logistic = CVtest( Group,
#                         function(x,...) glm(Group ~ ., mydat, , family=binomial(link="logit"), ...),
#                         function(obj, x) predict(obj, mydat, type='response') )
# 
# 
# resmultinom = CVtest(Group,
#                      function(x,...) multinom( Group ~ ., mydat , ...),
#                      function(obj, x) predict(obj, mydat, type="probs"), maxit=1000, trace=FALSE
# )

prob.lda = CVtest( Group, 
                   function(x,...) lda(Group ~ ., mydat, ...),
                   function(obj, x) predict(obj, mydat)$posterior )

res.logistic = ifelse(prob.logistic>0.5, 1,0)
# optimise threshold
# res.logistic = as.factor(ifelse(prob.logistic>OpThresh(prob.logistic, as.factor(Group), Grid=seq(0.25,0.75,by=0.02), Criteria="acc", positive="1"), 1, 0))
con(res.logistic, Group)
Acc = sum(res.logistic != Group)/length(Group)

top = opt10
layout(matrix(1:6, nrow=2))

my_colors <- adjustcolor(c(
  "#984ea3",
  "#377eb8",
  "#4daf4a", 
  "#a65628", 
  "#e41a1c"), alpha.f =0.5)

#Manuscript Figure 4B
for (ii in 1:length(top)) boxplot(mydat[,top[ii]] ~ pd$ConditionOrdered, main=markers[top[ii], "gene_symbol"], col = my_colors)
plot( roc(Group, prob.lda), print.auc=TRUE, print.thres=TRUE, main="LDA - LOOCV on train data" )

#plot( roc(Group, prob.logistic), print.auc=TRUE, print.thres=TRUE )
#plot( roc(Group, resmultinom), print.auc=TRUE, print.thres=TRUE )
#plot( roc(Group, prob.lda), print.auc=TRUE, print.thres=TRUE )
coords(roc(Group, prob.lda), "best", ret = c("threshold", "specificity", "sensitivity"))


##################################
# five-gene panel performance in test data
##################################

Group = as.numeric(as.factor(mydat.test$Status)) - 1
rand = 1:length(Group)
# mydat = mydat.test[,colnames(mydat.test) %in% opt10,drop=FALSE]
# prob.logistic = CVtest( Group,
#                         function(x,...) glm(Group ~ ., mydat, , family=binomial(link="logit"), ...),
#                         function(obj, x) predict(obj, mydat, type='response') )
# 
# 
# resmultinom = CVtest(Group,
#                      function(x,...) multinom( Group ~ ., mydat , ...),
#                      function(obj, x) predict(obj, mydat, type="probs"), maxit=1000, trace=FALSE
# )

prob.lda = CVtest( Group, 
                   function(x,...) lda(Group ~ ., mydat, ...),
                   function(obj, x) predict(obj, mydat)$posterior )

res.logistic = ifelse(prob.logistic>0.5, 1,0)
# optimise threshold
# res.logistic = as.factor(ifelse(prob.logistic>OpThresh(prob.logistic, as.factor(Group), Grid=seq(0.25,0.75,by=0.02), Criteria="acc", positive="1"), 1, 0))
con(res.logistic, Group)
Acc = sum(res.logistic != Group)/length(Group)

my_colors2 <- adjustcolor(c(
  "#984ea3",
  "#377eb8",
  "#4daf4a", 
  "#e41a1c"), alpha.f =0.5)

#Manuscript Figure 4C
top = opt10
layout(matrix(1:6, nrow=2))
for (ii in 1:length(top)) boxplot(mydat.test[,top[ii]] ~ test.metadata$ConditionOrdered, main=markers[top[ii], "gene_symbol"], col=my_colors2)
plot( roc(Group, prob.lda), print.auc=TRUE, print.thres=TRUE , main="LDA - LOOCV on test data")

coords(roc(Group, prob.lda), "best", ret = c("threshold", "specificity", "sensitivity"))
