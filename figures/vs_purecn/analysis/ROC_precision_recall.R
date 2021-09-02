library(data.table)
library(ggplot2)
library(MLmetrics)
library(ROCR)


auc_roc <- function(preds, actuals, returnDT=FALSE){
  # Calculate area under the ROC curve
  # If returnDT = TRUE, a data.table is returned
  
  #--------------------------------------------------
  # Hack to pass 'no visible binding for global variable' notes from R CMD check
  
  Pred <- NULL
  Actual <- NULL
  CumulativeFPR <- NULL
  CountFalse <- NULL
  CumulativeTPR <- NULL
  CountTrue <- NULL
  AdditionalArea <- NULL
  CumulativeArea <- NULL
  
  #--------------------------------------------------
  
  # Check if every prediction is identical and if so, return 0.5
  if(length(unique(preds)) == 1L) return(0.5)
  
  # Convert actuals to numeric if it's an ordered factor
  if(is(actuals, "factor")){
    if(is.ordered(actuals) & length(levels(actuals)) == 2) actuals <- as.numeric(actuals) - 1 else stop("actuals is type factor, but is unordered. Make it an ordered factor.")
  }
  
  dt <- data.table(Pred=preds, Actual=actuals*1L)
  setorder(dt, -Pred)
  
  dt <- dt[ , {
    CountTrue = sum(Actual)
    list(CountFalse=.N - CountTrue, CountTrue=CountTrue)
  }, by=Pred]
  
  # Calculate the CumulativeFalsePositiveRate and CumulativeTruePositiveRate
  dt[, CumulativeFPR := cumsum(CountFalse)/sum(CountFalse)]
  dt[, CumulativeTPR := cumsum(CountTrue)/sum(CountTrue)]
  
  # Calculate AUC ROC
  dt[, AdditionalArea := c(head(CumulativeFPR, 1) * head(CumulativeTPR, 1)/2,
                           (tail(CumulativeFPR, -1) - head(CumulativeFPR, -1)) * (head(CumulativeTPR, -1) + (tail(CumulativeTPR, -1) - head(CumulativeTPR, -1))/2))]
  dt[, CumulativeArea := cumsum(AdditionalArea)]
  
  # Return the desired result
  if(returnDT) return(dt[]) else return(tail(dt$CumulativeArea, 1))
}

ptm.tcga <- fread('../tcga/results/pure_tabnet_merge.csv')

data.splits <- fread('../../tmb/results/data_splits_210_with_tmb41.csv')


# for ROC or precision recall curve
get.tpr.fpr <- function(dt, percent, prob.column = "tabnet_proba_1"){
  ntp <- nrow(dt[get(prob.column) > percent][ truth == 1])
  nfp <- nrow(dt[get(prob.column) > percent][ truth == 0])
  nfn <- nrow(dt[get(prob.column) < percent][ truth == 1])
  ntn <- nrow(dt[get(prob.column) < percent][ truth == 0])
  # these are the three items needed for ROC plots and precision recall curves.
  # TPR (sensitivity aka recall) is mutual to both.  
  # FPR is ROC only.  Precision is P-R curve only.
  tpr <- ntp / (nfn + ntp) # sensitivity
  fpr <- nfp / (nfp + ntn) # fpr 
  precision <- ntp / (ntp + nfp) # precision 
  return(list(tpr,fpr, precision))
} 

get.tpr.fpr(ptm.tcga, 0.95, prob.column = 'tabnet_proba_1')
get.tpr.fpr(ptm.tcga, 0.95)
get.tpr.fpr(ptm.tcga, 0.95, prob.column = 'POSTERIOR.SOMATIC')


# takes a subset of the full venn merge table
build.curve.dt <- function(merged.dt, n_points = 100){
  out <- list()
  i <- 1
  tn.quantiles <- quantile(merged.dt$tabnet_proba_1, seq(0,1,1/n_points), na.rm = TRUE)
  pcn.quantiles <- quantile(merged.dt$POSTERIOR.SOMATIC, seq(0,1,1/n_points), na.rm = TRUE)
  for(point in 1:n_points){
    tn.prob.cutoff <- tn.quantiles[[point]]
    pcn.prob.cutoff <- pcn.quantiles[[point]]
    t.list <- get.tpr.fpr(merged.dt, tn.prob.cutoff, prob.column = 'tabnet_proba_1')
    pcn.list <- get.tpr.fpr(merged.dt, pcn.prob.cutoff, prob.column = 'POSTERIOR.SOMATIC')
    t.prob.cutoff.dt <- data.table(prob_cutoff = tn.prob.cutoff, tpr = t.list[[1]], fpr = t.list[[2]], precision = t.list[[3]], caller = 'tabnet' )
    pcn.prob.cutoff.dt <- data.table(prob_cutoff = pcn.prob.cutoff, tpr = pcn.list[[1]], fpr = pcn.list[[2]], precision = pcn.list[[3]], caller = 'purecn' )
    out[[i]] <- t.prob.cutoff.dt
    out[[i+ 1]] <- pcn.prob.cutoff.dt
    i <- i + 2
  }
  final <- rbindlist(out)
  final[,method := 'TabNet']
  final[caller == 'purecn',method := 'PureCN']
  final[,f1 := 2*tpr*precision/(tpr + precision)]
  return(final)
}

roc.dt <- build.curve.dt(ptm.tcga, 10)

roc.dt <- build.curve.dt(ptm.tcga, 100)

# see if it's working
ggplot(roc.dt, aes(x = fpr, y = tpr)) + geom_point(aes(color = caller))
ggplot(roc.dt, aes(x = tpr, y = precision)) + geom_point(aes(color = caller))

ptm.train <- ptm.tcga[sample %in% data.splits[split == 'Training set']$patient]
ptm.val <- ptm.tcga[sample %in% data.splits[split == 'Validation set']$patient]
ptm.test <- ptm.tcga[sample %in% data.splits[split == 'Test set']$patient]

ptm.hugo <- fread('../hugo/results/pure_tabnet_merge.csv')

theme.figure <- theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

roc.dt.train <- build.curve.dt(ptm.train, 500)
roc.dt.train.indels <- build.curve.dt(ptm.train[nchar(ref) != nchar(alt)], 500)
roc.dt.train.snvs <- build.curve.dt(ptm.train[nchar(ref) == nchar(alt)], 500)

roc.dt.val <- build.curve.dt(ptm.val, 500)

roc.dt.test <- build.curve.dt(ptm.test, 500)


# train PR curve
ggplot(roc.dt.train, aes(x = tpr, y = precision)) + geom_point(aes(color = method)) + xlab('recall') + ylim(0, 1) + theme.figure

# find optimal threshold for F1-score.
max(roc.dt.train[caller == 'purecn']$f1, na.rm = TRUE)
#0.8423194
roc.dt.train[caller == 'purecn'][f1 > 0.842] # roc.dt.train[caller == 'tabnet']$f1
#cutoff -- 0.005106713
# do 0.005 for cut-off in paper!!
max(roc.dt.train[caller == 'tabnet']$f1, na.rm = TRUE)
#0.8682862
roc.dt.train[caller == 'tabnet'][f1 > 0.868] # roc.dt.train[caller == 'tabnet']$f1
# optimal cutoff -- 0.4777132, close enough to 0.5

# SNVs and INDELs!
max(roc.dt.train.snvs[caller == 'tabnet']$f1, na.rm = TRUE)
# 0.8849
roc.dt.train.snvs[caller == 'tabnet'][f1 > 0.8849]
# 0.5075294
max(roc.dt.train.indels[caller == 'tabnet']$f1, na.rm = TRUE)
# 0.732... wowzers! :)
roc.dt.train.indels[caller == 'tabnet'][f1 > 0.732]
# 0.1368

# purecn snvs and indels
max(roc.dt.train.snvs[caller == 'purecn']$f1, na.rm = TRUE)
# 0.8445
roc.dt.train.snvs[caller == 'purecn'][f1 > 0.8445]
# 0.00495
max(roc.dt.train.indels[caller == 'purecn']$f1, na.rm = TRUE)
# 0.8119
roc.dt.train.indels[caller == 'purecn'][f1 > 0.8119]
# 0.014

max(roc.dt.val[caller == 'tabnet']$f1, na.rm = TRUE)
#0.6382356
roc.dt.val[caller == 'tabnet'][f1 > 0.638] # roc.dt.train[caller == 'tabnet']$f1
#cutoff - 0.548
max(roc.dt.val[caller == 'purecn']$f1, na.rm = TRUE)
# 0.6975
roc.dt.val[caller == 'purecn'][f1 > 0.697]
# cutoff - 0.01180693
# use optimal train cutoff -- 
roc.dt.val[caller == 'purecn'][prob_cutoff > 0.0053][prob_cutoff < 0.0057]
# 0.684, very close.


ggplot(roc.dt.test, aes(x = fpr, y = tpr)) + geom_point(aes(color = caller))
# :)
ggplot(roc.dt.test, aes(x = tpr, y = precision)) + geom_point(aes(color = method)) + xlab('recall') + ylim(0, 1) + theme.figure


max(roc.dt.test[caller == 'purecn']$f1, na.rm = TRUE)
min(roc.dt.test[caller == 'purecn']$f1, na.rm = TRUE)

max(roc.dt.test[caller == 'tabnet']$f1, na.rm = TRUE)
#0.885
max(roc.dt.test[caller == 'tabnet']$f1, na.rm = TRUE)

max(roc.dt.test[caller == 'purecn']$f1, na.rm = TRUE)



# load up time benchmark data
purecn.time <- fread('purecn_time_benchmark_results.csv')
tabnet.time <- fread('tabnet_time_benchmark_results.csv')
purecn.time[,method := 'PureCN']
purecn.time[,elapsed := purecn_elapsed]
tabnet.time[,method := 'TabNet']
tabnet.time[,elapsed := tabnet_elapsed]
benchmark.time <- rbind(purecn.time, tabnet.time, fill = TRUE)

ggplot(benchmark.time) +  geom_boxplot(aes(x = method, y  = elapsed, fill = method)) + theme.figure + ylab('Run time (s)')

wilcox.test(benchmark.time[method == 'PureCN']$elapsed, benchmark.time[method == 'TabNet']$elapsed)
mean.pcn.time <- mean(benchmark.time[method == 'PureCN']$elapsed)
mean.tn.time <- mean(benchmark.time[method == 'TabNet']$elapsed)

mean.pcn.time/mean.tn.time


# old stuff below



plot.roc <- function(dt){
    ggplot(dt, aes(x = fpr, y = tpr)) + geom_point(aes(color = caller)) + 
    xlab("False Positive Rate") + ylab("True Positive Rate") + theme_bw()
}

plot.pr <- function(dt, ymin = 0.4){
    ggplot(dt, aes(x = tpr, y = precision)) + geom_point(aes(color = caller)) +
    xlab("Recall") + ylab("Precision") + theme_bw() + ylim(ymin, 1.005)
}

plot.roc(roc.dt)

roc(ptm.ho$truth, ptm.ho$tabnet_proba_1)
# auc tabnet:  99.8
roc(ptm.ho$truth, ptm.ho$POSTERIOR.SOMATIC)
# auc purecn:  86.36


plot.pr(roc.dt)

# # with outliers removed 

ptm.no.outliers <- ptm.ho[sample %in% data.splits[actual.tmb <= 12]$V1]
length(unique(ptm.no.outliers$sample))

normal.tmb.roc.dt <- build.curve.dt(ptm.no.outliers, 100)

plot.pr(normal.tmb.roc.dt)
# save this

plot.roc(normal.tmb.roc.dt)
# save this


roc(ptm.no.outliers$truth, ptm.no.outliers$tabnet_proba_1)
# auc tabnet:  99.5
roc(ptm.no.outliers$truth, ptm.no.outliers$POSTERIOR.SOMATIC)
# auc purecn:  92.8


#get by indels

get.perf <- function(dt, percent, prob.column = "tabnet_proba_1", algo = "", dataset = "", v.cat = ""){
  ntp <- nrow(dt[get(prob.column) > percent][ truth == 1])
  nfp <- nrow(dt[get(prob.column) > percent][ truth == 0])
  nfn <- nrow(dt[get(prob.column) < percent][ truth == 1])
  ntn <- nrow(dt[get(prob.column) < percent][ truth == 0])
  sens <- ntp/(ntp + nfn)
  spec <- ntn/(ntn + nfp)
  ppv <- ntp/(ntp + nfp)
  npv <- ntn/(ntn + nfn)
  jaccard <- ntp/(nfp + ntp + nfn)
  f1 <- 2*(sens*ppv)/(sens+ppv)
  balanced.accuracy <- (sens + spec)/2
  total.called <- ntp + nfp + nfn + ntn
  call.rate <- total.called / nrow(dt[truth %in% c(0,1)])* 100
  # auc is agnostic to the supplied threshold
  auc <- auc_roc(get(prob.column, dt[!is.na(get(prob.column)) & !is.na(truth)]), dt[!is.na(get(prob.column)) & !is.na(truth)]$truth)
  ml.auc <- MLmetrics::AUC(get(prob.column, dt[!is.na(get(prob.column)) & !is.na(truth)]), dt[!is.na(get(prob.column)) & !is.na(truth)]$truth)
  ml.prauc <- MLmetrics::PRAUC(get(prob.column, dt[!is.na(get(prob.column)) & !is.na(truth)]), dt[!is.na(get(prob.column)) & !is.na(truth)]$truth)
  cat('sens',sens, 'spec',spec, 'ppv', ppv, 'npv', npv, 'jaccard', jaccard,'\n')
  cat('nFP', nfp, 'nTP', ntp, 'nFN',nfn, 'nTN',ntn,'\n')
  cat(sens, spec,  ppv,  npv,  jaccard,'\n')
  cat( nfp,  ntp, nfn, ntn)
  cat( auc)
  d1 <- data.table(Method = algo, "Dataset" = dataset, 'Variant category' = v.cat,
                   'AUC' = round(auc,3),
                   'mlmetrics.AUC' = round(ml.auc,3),
                   'mlmetrics.PRAUC' = round(ml.prauc,3),
                   Sensitivity = round(sens,3), Specificity = round(spec,3), 
                   'Positive predictive value' = round(ppv,3), 
                   'Negative predictive value' = round(npv,3), 
                   'F1-score' = round(f1,3), 
                   'Balanced accuracy' = round(balanced.accuracy,3),
                   'Jaccard index' = jaccard, 'FP' = nfp, 'TP' = ntp, 
                   'FN' = nfn, 'TN' = ntn, 'Total classified variants' = total.called, 
                   'Call rate' = round(call.rate, 1))
  return(d1)
}

# SNVS and INDELS
ptm.train.snvs <- ptm.train[(nchar(ref) == nchar(alt)) | ((nchar(REF)) == (nchar(ALT)))]
ptm.train.indels <- ptm.train[(nchar(ref) != nchar(alt)) | ((nchar(REF)) != (nchar(ALT)))]
ptm.val.snvs <- ptm.val[(nchar(ref) == nchar(alt)) | ((nchar(REF)) == (nchar(ALT)))]
ptm.val.indels <- ptm.val[(nchar(ref) != nchar(alt)) | ((nchar(REF)) != (nchar(ALT)))]
ptm.test.snvs <- ptm.test[(nchar(ref) == nchar(alt)) | ((nchar(REF)) == (nchar(ALT)))]
ptm.test.indels <- ptm.test[(nchar(ref) != nchar(alt)) | ((nchar(REF)) != (nchar(ALT)))]
ptm.hugo.snvs <- ptm.hugo[(nchar(ref) == nchar(alt)) | ((nchar(REF)) == (nchar(ALT)))]
ptm.hugo.indels <- ptm.hugo[(nchar(ref) != nchar(alt)) | ((nchar(REF)) != (nchar(ALT)))]

# TABNET, train, val, hugo
test.rocr <- function(dt, percent, prob.column = "tabnet_proba_1", algo = "", dataset = "", v.cat = ""){
  ml.auc <- MLmetrics::AUC(get(prob.column, dt[!is.na(get(prob.column)) & !is.na(truth)]), dt[!is.na(get(prob.column)) & !is.na(truth)]$truth)
  auc <- auc_roc(get(prob.column, dt[!is.na(get(prob.column)) & !is.na(truth)]), dt[!is.na(get(prob.column)) & !is.na(truth)]$truth)
  ml.prauc <- MLmetrics::PRAUC(get(prob.column, dt[!is.na(get(prob.column)) & !is.na(truth)]), dt[!is.na(get(prob.column)) & !is.na(truth)]$truth)
  return(ml.auc)
}

test.rocr(ptm.train, 0.5, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Training set', v.cat = 'Overall')

o.t.tr <- get.perf(ptm.train, 0.5, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Training set', v.cat = 'Overall')
s.t.tr <- get.perf(ptm.train.snvs, 0.508, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Training set', v.cat = 'SNVs')
i.t.tr <- get.perf(ptm.train.indels, 0.1368, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Training set', v.cat = 'Indels')
o.t.v <- get.perf(ptm.val, 0.5, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Validation set - COAD DLBC TGCT', v.cat = 'Overall')
s.t.v <- get.perf(ptm.val.snvs, 0.508, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Validation set - COAD DLBC TGCT', v.cat = 'SNVs')
i.t.v <- get.perf(ptm.val.indels, 0.1368, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Validation set - COAD DLBC TGCT', v.cat = 'Indels')
o.t.t <- get.perf(ptm.test, 0.5, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Blind test set - BRCA SARC UCEC', v.cat = 'Overall')
s.t.t <- get.perf(ptm.test.snvs, 0.508, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Blind test set - BRCA SARC UCEC', v.cat = 'SNVs')
i.t.t <- get.perf(ptm.test.indels, 0.1368, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Blind test set - BRCA SARC UCEC', v.cat = 'Indels')
o.t.h <- get.perf(ptm.hugo, 0.5, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Blind test set - metastatic melanoma', v.cat = 'Overall')
s.t.h <- get.perf(ptm.hugo.snvs, 0.508, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Blind test set - metastatic melanoma', v.cat = 'SNVs')
i.t.h <- get.perf(ptm.hugo.indels, 0.1368, prob.column = 'tabnet_proba_1', algo = 'TabNet', dataset = 'Blind test set - metastatic melanoma', v.cat = 'Indels')

# PURECN, train, val, hugo, using optiimal thresholds from training set.
o.p.tr <- get.perf(ptm.train, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Training set', v.cat = 'Overall')
s.p.tr <- get.perf(ptm.train.snvs, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Training set', v.cat = 'SNVs')
i.p.tr <- get.perf(ptm.train.indels, 0.014, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Training set', v.cat = 'Indels')
o.p.v <- get.perf(ptm.val, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Validation set - COAD DLBC TGCT', v.cat = 'Overall')
s.p.v <- get.perf(ptm.val.snvs, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Validation set - COAD DLBC TGCT', v.cat = 'SNVs')
i.p.v <- get.perf(ptm.val.indels, 0.014, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Validation set - COAD DLBC TGCT', v.cat = 'Indels')
o.p.t <- get.perf(ptm.test, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Blind test set - BRCA SARC UCEC', v.cat = 'Overall')
s.p.t <- get.perf(ptm.test.snvs, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Blind test set - BRCA SARC UCEC', v.cat = 'SNVs')
i.p.t <- get.perf(ptm.test.indels, 0.014, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Blind test set - BRCA SARC UCEC', v.cat = 'Indels')
o.p.h <- get.perf(ptm.hugo, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Blind test set - metastatic melanoma', v.cat = 'Overall')
s.p.h <- get.perf(ptm.hugo.snvs, 0.005, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Blind test set - metastatic melanoma', v.cat = 'SNVs')
i.p.h <- get.perf(ptm.hugo.indels, 0.014, prob.column = 'POSTERIOR.SOMATIC', algo = 'PureCN', dataset = 'Blind test set - metastatic melanoma', v.cat = 'Indels')

final.perf.dt <- rbindlist(list(o.t.tr,s.t.tr,i.t.tr,o.t.v, s.t.v, i.t.v, o.t.t, s.t.t, i.t.t, o.t.h, s.t.h, i.t.h, o.p.tr,s.p.tr, i.p.tr, o.p.v, s.p.v, i.p.v, o.p.t, s.p.t, i.p.t, o.p.h, s.p.h, i.p.h))
fwrite(final.perf.dt, file = 'final.accuracy.performance.csv')


