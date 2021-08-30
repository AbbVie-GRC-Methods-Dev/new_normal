library(ggplot2)
library(data.table)

data.splits <- data.table::fread('../../data_splits_210.csv')

full.c.maf <- data.table::fread('/unitynfs1/projects/mclaurt/TCGA_2020/tumor_only_caller/all_210_patients/tabnet/tabnet_pred_maf_popmax_filtered.csv')
#nrow: 165740
# variant level performance:
# output from tabnet
#split name: train
#TN: 48603, FN: 1401, TP : 20565, FP : 571
#sens: 0.936, spec: 0.988, ppv: 0.973, npv: 0.972
#Jaccard Index: 0.912
#split name: test
#TN: 23313, FN: 295, TP : 3761, FP : 362
#sens: 0.927, spec: 0.985, ppv: 0.912, npv: 0.988
#Jaccard Index: 0.851
#split name: validation
#TN: 35788, FN: 591, TP : 29865, FP : 625
#sens: 0.981, spec: 0.983, ppv: 0.980, npv: 0.984
#Jaccard Index: 0.961

full.c.maf[,result := '']
full.c.maf[(truth == 1) & (tabnet_pred == 1), result := 'TP']
full.c.maf[(truth == 0) & (tabnet_pred == 1), result := 'FP']
full.c.maf[(truth == 0) & (tabnet_pred == 0), result := 'TN']
full.c.maf[(truth == 1) & (tabnet_pred == 0), result := 'FN']

data.splits$split <- gsub("train", "Training set", data.splits$split)
# fix swapped validation / test nomenclature
data.splits$split <- gsub("test", "Validation set", data.splits$split)
data.splits$split <- gsub("validate", "Test set", data.splits$split)




four.metrics <- function(nFP, nTP, nFN, nTN ){
  sensitivity <- ifelse(((nTP + nFN ) > 0), nTP/(nTP + nFN),0)
  specificity <- ifelse(((nTN + nFP ) > 0), nTN/(nTN + nFP),0)
  ppv         <- ifelse(((nTP + nFP ) > 0),nTP/(nTP + nFP),0)
  npv         <- ifelse(((nTN + nFN ) > 0),nTN/(nTN + nFN),0)
  return(c(sensitivity, specificity, ppv, npv))
}

data.splits[, sens :=2] 
data.splits[, spec :=2]
data.splits[, ppv :=2]
data.splits[, npv :=2]
# without normal, with xgboost.
data.splits[, predicted.tmb :=-1]
# with normal.
data.splits[, actual.tmb :=-1]



process.maf.subset <- function(maf.subset, empty.data.splits){
   ds.copy <- copy(empty.data.splits)
   for(patient in unique(maf.subset$sample)){
     print(patient)
     patient.maf = maf.subset[sample == patient] 
     nTP = nrow(patient.maf[result == 'TP'])
     nFP = nrow(patient.maf[result == 'FP'])
     nFN = nrow(patient.maf[result == 'FN'])
     nTN = nrow(patient.maf[result == 'TN'])
     metrics <- four.metrics(nFP, nTP, nFN, nTN)           
     ds.copy[V1 == patient, sens := metrics[1]]
     ds.copy[V1 == patient, spec := metrics[2]]
     ds.copy[V1 == patient, ppv := metrics[3]]
     ds.copy[V1 == patient, npv := metrics[4]]
     ds.copy[V1 == patient, naive.tmb := ( nFP + nTP + nFN + nTN)/38]
     ds.copy[V1 == patient, predicted.tmb := (nFP + nTP)/38]
     ds.copy[V1 == patient, actual.tmb := (nTP + nFN)/38]
   }
   return(ds.copy)
}

perf.all <- process.maf.subset(full.c.maf, data.splits)
# inspect outliers on plots
#perf.all[sens < 0.1]
#perf.all[ppv < 0.2]
perf.all.normal.tmb <- perf.all[actual.tmb > 1][actual.tmb < 12] 

perf.indels <- process.maf.subset(full.c.maf[nchar(ref) != nchar(alt)], data.splits)
perf.snvs <- process.maf.subset(full.c.maf[nchar(ref) == nchar(alt)], data.splits)
# no negative values showing up in either of these!  so if the patient has no indels, for example,
# they will show up in the table and 
# have a sensitivity of 0 (16 patients are like this.)

#apply mean summary statistic to all patients.
# for averaging across patients
cols = sapply(perf.all, is.numeric)
cols = names(cols)[cols]
perf.all[, lapply(.SD, mean), .SDcols = cols]
perf.snvs[, lapply(.SD, mean), .SDcols = cols]

# indels only

# save this data.splits dt containing  tmb and TN/TP/FP/FN # to this folder.  
# for the Venn Diagrams with PureCN, this will be very useful.
fwrite(perf.all, file = 'data_splits_210_with_tmb38.csv')


for(patient in unique(data.splits$V1)){
   print(patient)
   patient.maf = full.c.maf[sample == patient] 
   nTP = nrow(patient.maf[result == 'TP'])
   nFP = nrow(patient.maf[result == 'FP'])
   nFN = nrow(patient.maf[result == 'FN'])
   nTN = nrow(patient.maf[result == 'TN'])
   metrics <- four.metrics(nFP, nTP, nFN, nTN)           
   data.splits[V1 == patient, sens := metrics[1]]
   data.splits[V1 == patient, spec := metrics[2]]
   data.splits[V1 == patient, ppv := metrics[3]]
   data.splits[V1 == patient, npv := metrics[4]]
   data.splits[V1 == patient, naive.tmb := ( nFP + nTP + nFN + nTN)/38]
   data.splits[V1 == patient, predicted.tmb := (nFP + nTP)/38]
   data.splits[V1 == patient, actual.tmb := (nTP + nFN)/38]
}
# reorder factors so they plot in the correct order.
data.splits$split_f = factor(data.splits$split, levels=c('Training set','Validation set','Test set'))
# handy function for showing the nonempty factors in each facet
boxplot_x_free <-facet_grid(~split_f, scales = "free_x",space = 'free')

# make black and white and remove the annoying legend (it's redundant in these boxplots)
tcga_boxplot_theme <- theme_bw() + theme(legend.position = 'none',strip.background = element_rect(fill = 'white'))

boxplot.metric.tissue.split <- function(metric,label, log.y = FALSE){
  p <- ggplot(data.splits) + geom_boxplot(aes(x = tissue, y = get(metric), fill = tissue)) + ylab(label) + boxplot_x_free + tcga_boxplot_theme 
  if(log.y){p <- p + scale_y_continuous(trans='log2')}
  p
}

# sensitivity, specificity results broken down by indication.
const.h = 2.5
boxplot.metric.tissue.split('sens',label = 'Sensitivity')
ggsave('sensitivity.indication.boxplot.png', height = const.h, width = const.h*2.5, units = 'in')
boxplot.metric.tissue.split('spec', label = 'Specificity')
ggsave('specificity.indication.boxplot.png', height = const.h, width = const.h*2.5, units = 'in')
boxplot.metric.tissue.split('ppv', label = 'Positive predictive value')
ggsave('ppv.indication.boxplot.png', height = const.h, width = const.h*2.5, units = 'in')
boxplot.metric.tissue.split('npv', label = 'Negative predictive value')
ggsave('npv.indication.boxplot.png', height = const.h, width = const.h*2.5, units = 'in')
boxplot.metric.tissue.split('actual.tmb', label = 'Mutational burden (with matched normal)', log.y = TRUE)
ggsave('actualtmb.indication.boxplot.png', height = const.h, width = const.h*2.5, units = 'in')

# simple tmb plots.  no R squared or lm fit.  Just exploratory because i need to determine an aactual.tmb threshold for plotting.

simple.3fold.scatter.plot <- function(actual.tmb.upper.bound = 10000, mode = 'predicted', log.axes = FALSE){
  if(mode == 'predicted'){
    p <- ggplot(data.splits[actual.tmb < actual.tmb.upper.bound]) + geom_point(aes(x = predicted.tmb, y = actual.tmb, color = tissue)) + theme_bw() + facet_grid(~split_f)
  }else if(mode == 'naive'){
    p <- ggplot(data.splits[actual.tmb < actual.tmb.upper.bound]) + geom_point(aes(x = naive.tmb, y = actual.tmb, color = tissue)) + theme_bw() + facet_grid(~split_f)
  }
  if(log.axes){p <- p + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')}
  p
}

# tumor-only predicted TMB vs actual TMB plots
simple.3fold.scatter.plot()
simple.3fold.scatter.plot(32)
simple.3fold.scatter.plot(12)

# tumor-only naive TMB vs actual plots
simple.3fold.scatter.plot( mode = 'naive')
simple.3fold.scatter.plot(32, mode = 'naive')
simple.3fold.scatter.plot(12, mode = 'naive')

get.lm.summary <- function(mode = 'predicted', dt){
  if(mode == 'predicted'){
    return(summary(lm(actual.tmb ~ predicted.tmb, data = dt)))
  }else if(mode == 'naive'){
    return(summary(lm(actual.tmb ~ naive.tmb, data = dt)))
  }
}

# overall for all three, then just on test / validation for ML and naive.

overall.s <- get.lm.summary(dt = data.splits)
predicted.s <- get.lm.summary(dt = data.splits[split_f %in% c('Validation set','Test set')])
test.val.naive.s <- get.lm.summary(mode = 'naive', dt = data.splits[split_f %in% c('Validation set','Test set')])


get.pval <- function(lm.model.summary){
  fstat <- lm.model.summary$fstat
  p.val <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  return(p.val)
}


train.test.val.fit.dt <- function(thresholded.dt){
  fit.dt <- data.table() 
  for(s in c('Training set', 'Validation set', 'Test set')){
    # get lm fit for tumor only prediction of TMB 
    pred.summary <- get.lm.summary(dt = thresholded.dt[split_f == s])
    ttv.pred.dt <- data.table(mode = 'predicted', split_f = s, slope = pred.summary$coefficients['predicted.tmb','Estimate'], r.sq = pred.summary$r.squared, p.val = get.pval(pred.summary))
    # get lm fit for naive TMB estimation
    naive.summary <- get.lm.summary(mode = 'naive', dt = thresholded.dt[split_f == s])
    ttv.naive.dt <- data.table(mode = 'naive',split_f = s, slope = naive.summary$coefficients['naive.tmb','Estimate'], r.sq = naive.summary$r.squared, p.val = get.pval(naive.summary))
    # should usually not be rbinding in a for loop, but this is a very
    # small number of loops!
    fit.dt <- rbind(fit.dt, ttv.pred.dt, ttv.naive.dt)
  }
  
  fit.dt[, fit.label := paste0('slope: ', round(slope, digits = 3),'\nR-squared: ',round(r.sq, digits = 3), ',\np-value: ', formatC(p.val, format = "e", digits = 2))]
  # get levels in correct order.
  fit.dt$split_f = factor(fit.dt$split_f, levels=c('Training set','Validation set','Test set'))
  return(fit.dt)
}

# all TMBs.
# n = 210.
# :)
ttv.fit.dt <- train.test.val.fit.dt(data.splits)

# fit for samples with true TMB <= 12.
# n = 193.
# :) :) :)
ttv.fit.dt.lte12 <- train.test.val.fit.dt(data.splits[actual.tmb <= 12])





pred.fit.dt <- data.table(split_f =  'all predicted',  r.sq= predicted.s$r.squared, p.val = get.pval(predicted.s))
pred.fit.dt[,fit.label := paste0('R-squared: ',round(r.sq, digits = 3), ',\np-value: ', formatC(p.val, format = "e", digits = 2))]


max(data.splits$actual.tmb)
#[1] 327.2368

# tmb threshold?
tmb.ecdf <- ecdf(data.splits$actual.tmb)
plot(tmb.ecdf)
tmb.ecdf <- ecdf(data.splits[actual.tmb < 100]$actual.tmb)
plot(tmb.ecdf)
tmb.ecdf <- ecdf(data.splits[actual.tmb < 20]$actual.tmb)
plot(tmb.ecdf)
# looks like 5.
table(data.splits[actual.tmb > 12]$tissue)
#BLCA COAD HNSC LUAD LUSC STAD UCEC 
#1    1    1    5    1    1    7 

tmb.triptych <- function(actual.tmb.threshold = max(data.splits$actual.tmb), mode = 'predicted', log.axes = FALSE, scales.option ='fixed'){
  mode.name <- mode
  input.data <- data.splits[actual.tmb <= actual.tmb.threshold]
  lm.fit.dt <- train.test.val.fit.dt(input.data)[mode == mode.name]
  print(lm.fit.dt) 
  x.axis <- ifelse(mode.name == 'predicted', 'predicted.tmb','naive.tmb')
  x.lab <- ifelse(mode.name == 'predicted', 'Tumor-Only Predicted TMB','Tumor-Only Naive TMB')
  
  max.x.axis.tmb <- max(get(x.axis,input.data))
  p <- ggplot(input.data, aes(x = get(x.axis), y = actual.tmb )) + geom_point( aes(color = tissue))
  p <- p + geom_smooth(method = 'lm', color = 'black') + theme_bw() + theme(strip.background = element_rect(fill = 'white')) +  facet_grid(~split_f, scales = scales.option)
  p <- p + ylim(0, actual.tmb.threshold) + xlab(x.lab) + ylab('Matched Normal TMB') 
  if(scales.option == 'fixed'){
    p <- p# + xlim(0, max.x.axis.tmb)
  }
  if(log.axes){
      p <- p + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')
      p <- p + geom_text(data = lm.fit.dt, aes(x = max.x.axis.tmb*.035, y = actual.tmb.threshold*.85, label = fit.label) )
  }else{
    if(scales.option == 'fixed'){
      p <- p + geom_text(data = lm.fit.dt, aes(x = max.x.axis.tmb *.35, y = actual.tmb.threshold*.85, label = fit.label) )
      } else{
      p <- p + geom_text(data = lm.fit.dt, aes(x = 22.35, y = actual.tmb.threshold*.85, label = fit.label) )
    }
    }
  p
}

const.h2 <- 4.5

tmb.triptych(mode = 'naive')
ggsave('tmb.naive.lin.png', height = const.h2, width = const.h2*2, units = 'in')
tmb.triptych(mode = 'naive', log.axes = TRUE)
ggsave('tmb.naive.log.png', height = const.h2, width = const.h2*2, units = 'in')

tmb.triptych(actual.tmb.threshold = 12, mode = 'naive', scales.option = 'free_x')
#tmb.triptych(actual.tmb.threshold = 12, mode = 'naive')
ggsave('tmb.naive.lin12.png', height = const.h2, width = const.h2*2, units = 'in')

tmb.triptych(mode = 'predicted')
ggsave('tmb.predicted.lin.png', height = const.h2, width = const.h2*2, units = 'in')
tmb.triptych(mode = 'predicted', log.axes = TRUE)
ggsave('tmb.predicted.log.png', height = const.h2, width = const.h2*2, units = 'in')
tmb.triptych(actual.tmb.threshold = 12, mode = 'predicted')
ggsave('tmb.predicted.lin12.png', height = const.h2, width = const.h2*2, units = 'in')
tmb.triptych(actual.tmb.threshold = 12, mode = 'predicted', log.axes = TRUE)
ggsave('tmb.predicted.log12.png', height = const.h2, width = const.h2*2, units = 'in')

# ok, full linear, TMB-12 cutoff, and log log plots all saved.  :)
