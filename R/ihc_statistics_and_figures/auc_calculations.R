library(pROC)

# AUC calculation for IHC

dat = read.csv("all_scores_edited.csv")

dd = dat
pd = data.frame(Sample=colnames(dat)[-c(1)])  
pd$Condition = substring(pd$Sample, 1,1)
table(pd$Condition)

pd$Status = 1
pd$Status[pd$Condition %in% c("N")] = 0

mydat = data.frame(as.numeric(as.factor(pd$Status))-1, t(dd[,-1]))
colnames(mydat) = c("Status", paste0(dd[,1]))

gp = mydat$Status
vec = mydat$X

roc.res = roc(gp, vec)
# change direction for MLN ">", if directions is defines, default = auto
names(roc.res)  # includes AUC
plot(roc.res, print.auc=TRUE, main="X", asp=NA, print.thres="best", print.thres.best.method="youden")
ci.auc(roc.res)

auc_ci <- ci.auc(roc.res)

plot(roc.res, print.auc=FALSE, main="X", asp=NA, print.thres=TRUE)
text(
  0.4, 0.2,
  labels = paste0(
    "AUC = ", round(auc(roc.res), 3),
    "\n95% CI: ",
    round(auc_ci[1], 3), "-",
    round(auc_ci[3], 3)
  )
)

coords(
  roc.res,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")
)
