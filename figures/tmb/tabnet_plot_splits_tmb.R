library(ggplot2)
library(data.table)

tcga.dir <- '../../test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/'
data.splits <- data.table::fread(paste0(tcga.dir,'data_splits.csv'), header = TRUE)

full.c.maf <- data.table::fread(paste0(tcga.dir, 'mafs/tabnet_trained_tcga_v3_predictions_all_features.csv'))
 
full.c.maf[,result := '']
full.c.maf[(truth == 1) & (tabnet_pred == 1), result := 'TP']
full.c.maf[(truth == 0) & (tabnet_pred == 1), result := 'FP']
full.c.maf[(truth == 0) & (tabnet_pred == 0), result := 'TN']
full.c.maf[(truth == 1) & (tabnet_pred == 0), result := 'FN']

data.splits$split <- gsub("train", "Training set", data.splits$split)
data.splits$split <- gsub("validation", "Validation set", data.splits$split)
data.splits$split <- gsub("test", "Test set", data.splits$split)


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
data.splits[, f1 :=2]
# without normal, with tabnet.
data.splits[, predicted.tmb :=-1]
# with normal.
data.splits[, actual.tmb :=-1]


process.maf.subset <- function(maf.subset, empty.data.splits){
   ds.copy <- copy(empty.data.splits)
   for(p in unique(maf.subset$sample)){
     print(p)
     patient.maf = maf.subset[sample == p] 
     nTP = nrow(patient.maf[result == 'TP'])
     nFP = nrow(patient.maf[result == 'FP'])
     nFN = nrow(patient.maf[result == 'FN'])
     nTN = nrow(patient.maf[result == 'TN'])
     # calculate median t_alt_freq for all true somatic mutations for this patient. 
     # akin to purity.
     p.median.true.somatic.t.alt.freq <- median(patient.maf[result %in% c('TP','FN')]$t_alt_freq)
     metrics <- four.metrics(nFP, nTP, nFN, nTN)           
     ds.copy[patient == p, sens := metrics[1]]
     ds.copy[patient == p, spec := metrics[2]]
     ds.copy[patient == p, ppv := metrics[3]]
     ds.copy[patient == p, npv := metrics[4]]
     ds.copy[patient == p, f1 := 2*sens*ppv/(sens + ppv)]
     # using 41, the weighted average of capture kit sizes.
     ds.copy[patient == p, naive.tmb := ( nFP + nTP + nFN + nTN)/41]
     ds.copy[patient == p, predicted.tmb := (nFP + nTP)/41]
     ds.copy[patient == p, actual.tmb := (nTP + nFN)/41]
     ds.copy[patient == p, median.true.somatic.t.alt.freq := p.median.true.somatic.t.alt.freq]
     
     
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
perf.indels[, lapply(.SD, mean), .SDcols = cols]
perf.all[split == 'Test set', lapply(.SD, mean), .SDcols = cols]

perf.all.normal.tmb[split == 'Test set', lapply(.SD, mean), .SDcols = cols]
perf.all.normal.tmb[split == 'Training set', lapply(.SD, mean), .SDcols = cols]
perf.all.normal.tmb[split == 'Validation set', lapply(.SD, mean), .SDcols = cols]

# indels only

# save this data.splits dt containing  tmb and TN/TP/FP/FN # to this folder.  
# for the Venn Diagrams with PureCN, this will be very useful.
fwrite(perf.all, file = 'data_splits_210_with_tmb41.csv')
fwrite(perf.snvs, file = 'snvs_data_splits_210_with_tmb41.csv')
fwrite(perf.indels, file = 'indels_data_splits_210_with_tmb41.csv')


# reorder factors so they plot in the correct order.
perf.all$split_f = factor(perf.all$split, levels=c('Training set','Validation set','Test set'))
# handy function for showing the nonempty factors in each facet
boxplot_x_free <-facet_grid(~split_f, scales = "free_x",space = 'free')

# make black and white and remove the annoying legend (it's redundant in these boxplots)
tcga_boxplot_theme <- theme_bw() + theme(legend.position = 'none',strip.background = element_rect(fill = 'white'))

boxplot.metric.tissue.split <- function(metric, label, log.y = FALSE){
  p <- ggplot(perf.all) + geom_boxplot(aes(x = indication, y = get(metric), fill = indication)) + xlab('Tissue') + ylab(label) + boxplot_x_free + tcga_boxplot_theme 
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
boxplot.metric.tissue.split('f1',label = 'F1-score')

# hmm, facet by metric?

perf.all.m <- melt(perf.all, measure.vars = c('sens','spec','ppv','npv'), variable.name = 'metric' )
ggplot(perf.all.m) + geom_boxplot(aes(x = indication, y = value, fill = indication)) + xlab('Tissue') + boxplot_x_free + facet_grid(metric~split_f, scales = "free_x",space = 'free') + tcga_boxplot_theme   
# better colors
# four beach colors
four.beach.colors <- c("#6C854D","#70C2C6", "#D2D9E8", "#cfc1a9")

ggplot(perf.all.m) + geom_boxplot(aes(x = indication, y = value, fill = metric)) + xlab('Tissue') + boxplot_x_free + facet_grid(~split_f, scales = "free_x",space = 'free') + tcga_boxplot_theme + scale_fill_manual(values = four.beach.colors)
ggsave('big.plot.pdf', height = 4, width = 7)

# might want to refine the training... looks like a bit of overfitting. 
# maybe harder indications?
ggplot(perf.all.m) + geom_boxplot(aes(x = indication, y = value, fill = metric)) + xlab('Tissue') + boxplot_x_free + facet_grid(~split_f, scales = "free_x",space = 'free') + tcga_boxplot_theme + scale_fill_manual(values = c("#6C854D","#70C2C6", "#D2D9E8", "#cfc1a9"))
ggsave('big.plot.pdf', height = 4, width = 7)

perf.all.m <- melt(perf.all[actual.tmb > 1][actual.tmb < 12] , measure.vars = c('sens','spec','ppv','npv'), variable.name = 'metric' )
ggplot(perf.all.m) + geom_boxplot(aes(x = indication, y = value, fill = metric)) + xlab('Tissue') + boxplot_x_free + facet_grid(~split_f, scales = "free_x",space = 'free') + tcga_boxplot_theme   
perf.all.m <- melt(perf.all[actual.tmb > 1] , measure.vars = c('sens','spec','ppv','npv'), variable.name = 'metric' )

ggplot(perf.all.m) + geom_boxplot(aes(x = indication, y = value, fill = metric)) + xlab('Tissue') + boxplot_x_free + facet_grid(~split_f, scales = "free_x",space = 'free') + tcga_boxplot_theme   

### look at hugo dataset file..
# ./hugo_performance.R

p <- rbind(perf.all,perf.h, fill = TRUE)
# plot sensitivity vs median somatic t.alt.count for example
# scatter plot
ggplot(perf.all) + geom_point(aes(x = median.true.somatic.t.alt.freq, y = sens, color = indication)) + theme_bw()
# with hugo and linear fit :)
ggplot(p, aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + 
      geom_smooth(method = 'lm', color = '#000000', se = FALSE ) + theme_bw() + 
      ylab('Model Sensitivity (TPR)') + xlab('Median True Somatic Variant Allele Fraction')

ggsave('sensitivity_vs_true_somatic_vaf.pdf', height = 4, width = 5.5)

ggplot(p[median.true.somatic.t.alt.freq <= 0.5], aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + geom_smooth(method = 'lm', color = '#000000', se = FALSE ) + theme_bw()
ggplot(p[indication == 'TGCT'], aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + geom_smooth(method = 'lm', color = '#000000', se = FALSE ) + theme_bw()
ggplot(p[indication == 'melanoma'], aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + geom_smooth(method = 'lm', color = '#000000', se = FALSE ) + ylab('TPR') + xlab('Median Somatic VAF') +  theme_bw()
ggplot(p[indication == 'GBM'], aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + geom_smooth(method = 'lm', color = '#000000', se = FALSE ) + ylab('TPR') + xlab('Median Somatic VAF') + theme_bw()
ggplot(p[indication == 'BLCA'], aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + geom_smooth(method = 'lm', color = '#000000', se = FALSE ) + ylab('TPR') +  xlab('Median Somatic VAF') + theme_bw()
# useful!  maybe as supplement. melanoma and GBM  among the most negative covariances.

t.alt.freq.sens.fit <- lm(sens ~ median.true.somatic.t.alt.freq, data = p)
summary(t.alt.freq.sens.fit)
# R-squared = 0.4
# quadratic
#t.alt.freq.sens.fit <- lm(sens ~ median.true.somatic.t.alt.freq + median.true.somatic.t.alt.freq^2, data = p)
#summary(t.alt.freq.sens.fit)
t.alt.freq.sens.fit <- lm(sens ~ median.true.somatic.t.alt.freq, data = p[median.true.somatic.t.alt.freq < 0.5])
summary(t.alt.freq.sens.fit)

# covariance and correlation of sensitivity dependence by tissue
library(boot)

cov.dt <- data.table()
boot.cov.dt <- data.table()
for(t in unique(p$indication)){
  t.subset.dt <- p[indication == t]
  t.cov <- cov(t.subset.dt$sens, t.subset.dt$median.true.somatic.t.alt.freq)
  t.cor <- cor(t.subset.dt$sens, t.subset.dt$median.true.somatic.t.alt.freq)
  cov.dt <- rbind(cov.dt, data.table(covariance = t.cov, correlation = t.cor, tissue = t))
  booted.cor <- boot(as.data.frame(t.subset.dt), 
             statistic = function(data, i) {
               cor(data[i, "median.true.somatic.t.alt.freq"], data[i, "sens"], method='pearson')
             },
             R = 1000
  )
  booted.cov <- boot(as.data.frame(t.subset.dt), 
             statistic = function(data, i) {
               cov(data[i, "median.true.somatic.t.alt.freq"], data[i, "sens"], method='pearson')
             },
             R = 1000
  )
  boot.cov.dt <- rbind(boot.cov.dt, data.table(tissue = t, boot.cor = booted.cor$t, boot.cov = booted.cov$t))
}

ggplot(cov.dt) + geom_bar(aes(x = tissue, y = correlation), stat = 'identity')
boot.cov.dt[,median.cov := median(boot.cov.V1) , .(tissue)]
boot.cov.dt[,median.cor := median(boot.cor.V1), .(tissue)]

ggplot(boot.cov.dt) + geom_boxplot(aes(x = reorder(tissue, boot.cov.V1, FUN = median), y = boot.cov.V1, fill = median.cov )) + xlab('Tissue') + ylab('Covariance between TPR and median somatic VAF') +scale_fill_viridis_c(alpha = 0.2)
ggplot(boot.cov.dt) + geom_boxplot(aes(x = reorder(tissue, boot.cor.V1, FUN = median), y = boot.cor.V1, fill = median.cor)) + xlab('Tissue') + ylab('Correlation between TPR and median somatic VAF') + scale_fill_viridis_c(alpha = 0.2)

ggplot(p, aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + geom_smooth(method = 'lm') + theme_bw()

ggplot(p, aes(x = median.true.somatic.t.alt.freq, y = sens)) + geom_point(aes(color = indication)) + geom_smooth(method = 'lm') + theme_bw()

# the higher the log-TMB, the higher the PPV.
ggplot(perf.all) + geom_point(aes(x = log(actual.tmb), y = ppv, color = indication)) + ylab('PPV') + theme_bw() + xlab("log (MN-TMB)")
# plus hugo melanoma
ggplot(p) + geom_point(aes(x = log(actual.tmb), y = ppv, color = indication)) + ylab('PPV') + theme_bw() + xlab("log (MN-TMB)")
ggsave('ppv_vs_log_mn-tmb.pdf', height = 4, width = 5.5)

# plus line of best fit.
ggplot(p,aes(x = log(actual.tmb), y = ppv)) + geom_point(aes(color = indication)) + ylab('PPV') + geom_smooth(method = 'lm', color = 'black', se = FALSE) + theme_bw() + xlab("log (MN-TMB)") + ylim(0, 1)

ppv.logtmb.fit <- lm(ppv ~ log(actual.tmb), data = p)

summary(ppv.logtmb.fit)

t_alt_freq_hist <- function(maf){
  ggplot(maf) + geom_histogram(aes(x = t_alt_freq, fill = factor(truth) )) + facet_wrap(~ sample) +  scale_fill_manual(values = c("#DD7089", "#3F92F6")) + theme_bw()
} 

t_alt_freq_hist(full.c.maf[sample %in% perf.all$patient[1:20] ])


# LAML
t_alt_freq_hist(full.c.maf[sample %in% perf.all[indication == 'LAML']$patient])

# BRCA
t_alt_freq_hist(full.c.maf[sample %in% perf.all[indication == 'BRCA']$patient])

#TGCT
t_alt_freq_hist(full.c.maf[sample %in% perf.all[indication == 'TGCT']$patient])

#UCEC
t_alt_freq_hist(full.c.maf[sample %in% perf.all[indication == 'UCEC']$patient])

#BRCA
t_alt_freq_hist(full.c.maf[sample %in% perf.all[indication == 'BRCA']$patient])

# negative results
#library(diptest)
#hist(full.c.maf[,dip(t_alt_freq), sample]$V1)
#dips <- full.c.maf[,dip(t_alt_freq), sample]
#setnames(old = 'sample', new = 'patient', dips)
#setnames(old = 'V1', new = 'dip_score', dips)
#m  <- merge(perf.all.no.laml, dips, by ="patient")
#ggplot(m) + geom_point(aes(x = dip_score, y = sens))
#ggplot(m) + geom_point(aes(x = dip_score, y = f1))
# no relationship



# simple tmb plots.  no R squared or lm fit.  Just exploratory because i need to determine an aactual.tmb threshold for plotting.

simple.3fold.scatter.plot <- function(perf.df, actual.tmb.upper.bound = 10000, mode = 'predicted', log.axes = FALSE){
  if(mode == 'predicted'){
    p <- ggplot(perf.df[actual.tmb < actual.tmb.upper.bound]) + geom_point(aes(x = predicted.tmb, y = actual.tmb, color = indication)) + theme_bw() + facet_grid(~split_f)
  }else if(mode == 'naive'){
    p <- ggplot(perf.df[actual.tmb < actual.tmb.upper.bound]) + geom_point(aes(x = naive.tmb, y = actual.tmb, color = indication)) + theme_bw() + facet_grid(~split_f)
  }
  if(log.axes){p <- p + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')}
  p
}

# tumor-only predicted TMB vs actual TMB plots
simple.3fold.scatter.plot(perf.all)
simple.3fold.scatter.plot(perf.all,32)
simple.3fold.scatter.plot(perf.all,12)

# tumor-only naive TMB vs actual plots
simple.3fold.scatter.plot(perf.all, mode = 'naive')
simple.3fold.scatter.plot(perf.all, 32, mode = 'naive')
simple.3fold.scatter.plot(perf.all, 12, mode = 'naive')

get.lm.summary <- function(mode = 'predicted', dt){
  if(mode == 'predicted'){
    return(summary(lm(actual.tmb ~ predicted.tmb, data = dt)))
  }else if(mode == 'naive'){
    return(summary(lm(actual.tmb ~ naive.tmb, data = dt)))
  }
}

# overall for all three, then just on test / validation for ML and naive.

overall.s <- get.lm.summary(dt = perf.all)
predicted.s <- get.lm.summary(dt = perf.all[split_f %in% c('Validation set','Test set')])
test.val.naive.s <- get.lm.summary(mode = 'naive', dt = perf.all[split_f %in% c('Validation set','Test set')])


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
ttv.fit.dt <- train.test.val.fit.dt(perf.all)

# fit for samples with true TMB <= 12.
# n = 193.
# :) :) :)
ttv.fit.dt.lte12 <- train.test.val.fit.dt(perf.all[actual.tmb <= 12])





pred.fit.dt <- data.table(split_f =  'all predicted',  r.sq= predicted.s$r.squared, p.val = get.pval(predicted.s))
pred.fit.dt[,fit.label := paste0('R-squared: ',round(r.sq, digits = 3), ',\np-value: ', formatC(p.val, format = "e", digits = 2))]


max(perf.all$actual.tmb)
#[1] 300.39 # old:  324.1053

# tmb threshold?
tmb.ecdf <- ecdf(perf.all$actual.tmb)
plot(tmb.ecdf)
tmb.ecdf <- ecdf(perf.all[actual.tmb < 100]$actual.tmb)
plot(tmb.ecdf)
tmb.ecdf <- ecdf(perf.all[actual.tmb < 20]$actual.tmb)
plot(tmb.ecdf)
# looks like 5.
table(perf.all[actual.tmb > 12]$indication)
#BLCA COAD HNSC LUAD LUSC STAD UCEC 
#1    1    1    5    1    1    7 

tmb.triptych <- function(actual.tmb.threshold = max(perf.all$actual.tmb), mode = 'predicted',
                         log.axes = FALSE, scales.option ='fixed', labels.on = TRUE){
  mode.name <- mode
  input.data <- perf.all[actual.tmb <= actual.tmb.threshold]
  lm.fit.dt <- train.test.val.fit.dt(input.data)[mode == mode.name]
  print(lm.fit.dt) 
  x.axis <- ifelse(mode.name == 'predicted', 'predicted.tmb','naive.tmb')
  x.lab <- ifelse(mode.name == 'predicted', 'Tumor-Only Predicted TMB','Tumor-Only Naive TMB')
  
  max.x.axis.tmb <- max(get(x.axis,input.data))
  p <- ggplot(input.data, aes(x = get(x.axis), y = actual.tmb )) + geom_point( aes(color = indication))
  p <- p + geom_smooth(method = 'lm', color = 'black', se = FALSE) + theme_bw() + theme(strip.background = element_rect(fill = 'white')) +  facet_grid(~split_f, scales = scales.option)
  p <- p + ylim(0, actual.tmb.threshold) + xlab(x.lab) + ylab('Matched Normal TMB') 
  if(scales.option == 'fixed'){
    p <- p# + xlim(0, max.x.axis.tmb)
  }
  if(log.axes){
      p <- p + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')
      if(labels.on){
      p <- p + geom_text(data = lm.fit.dt, aes(x = max.x.axis.tmb*.035, y = actual.tmb.threshold*.85, label = fit.label) )
      }
  }else if (labels.on){
    if(scales.option == 'fixed'){
      p <- p + geom_text(data = lm.fit.dt, aes(x = max.x.axis.tmb *.15, y = actual.tmb.threshold*.85, label = fit.label) )
      } else{
      p <- p + geom_text(data = lm.fit.dt, aes(x = 22.35, y = actual.tmb.threshold*.85, label = fit.label) )
      }
  }
  #p <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(panel.border=element_blank(),
             strip.text=element_text(size=12, colour="black"),
             axis.line = element_line(colour = 'black'),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             strip.background=element_rect(colour="white", 
                                           fill="white"))
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

# :)
tmb.triptych(actual.tmb.threshold = 12, mode = 'naive', scales.option = 'free_x', labels.on = FALSE)
ggsave('tmb.naive.lin12.no.labels.pdf', height = const.h2, width = const.h2*2, units = 'in')

tmb.triptych(mode = 'predicted')
ggsave('tmb.predicted.lin.png', height = const.h2, width = const.h2*2, units = 'in')
tmb.triptych(mode = 'predicted', log.axes = TRUE)
ggsave('tmb.predicted.log.png', height = const.h2, width = const.h2*2, units = 'in')

tmb.triptych(actual.tmb.threshold = 12, mode = 'predicted')
ggsave('tmb.predicted.lin12.png', height = const.h2, width = const.h2*2, units = 'in')

# :)
tmb.triptych(actual.tmb.threshold = 12, mode = 'predicted', labels.on = FALSE)
ggsave('tmb.predicted.lin12.no.labels.pdf', height = const.h2, width = const.h2*2, units = 'in')

tmb.triptych(actual.tmb.threshold = 12, mode = 'predicted', log.axes = TRUE)
ggsave('tmb.predicted.log12.png', height = const.h2, width = const.h2*2, units = 'in')

# ok, full linear, TMB-12 cutoff, and log log plots all saved.  :)

# load up TCGA clinical data, association with naive-TMB and race?
tcga.clin <- data.table::fread('/Users/mclaurt/unitynfs1/projects/masicdx/tcga/clinical/clinical')

# relevant columns:
#"jewish_origin"   "ethnicity"   "race" "bcr_patient_barcode"

# this includes the training set
#merged.with.race <- merge(perf.all, tcga.clin[,.(bcr_patient_barcode, race, ethnicity, jewish_origin)], by.x = 'patient', by.y = 'bcr_patient_barcode')
# this does not include training set (or Hugo dataset.)
merged.with.race <- merge(p[split_f != 'Training set'], tcga.clin[,.(bcr_patient_barcode, race, ethnicity, jewish_origin)], by.x = 'patient', by.y = 'bcr_patient_barcode')
merged.with.race[, race.f := 'Unknown']
merged.with.race[race == "WHITE" , race.f := 'White']
merged.with.race[race == "BLACK OR AFRICAN AMERICAN" , race.f := 'Black']

table(merged.with.race$race)

b.vs.w.plot <- function(dt, mode){
  p <- ggplot(dt[race.f %in% c("Black", "White")], aes(x = race.f, fill = race.f))  + xlab('Race')
  if(mode == 'naive'){
    p <- p + geom_boxplot(aes(y = naive.tmb)) + ylab('Naive Tumor-Only TMB') 
  }else if(mode == 'actual'){
    p <- p + geom_boxplot(aes(y = actual.tmb)) + ylab('Actual Matched-Normal TMB') 
  }else if(mode == 'predicted'){
    p <- p + geom_boxplot(aes(y = predicted.tmb)) + ylab('ML-Predicted Tumor-Only TMB') 
  }
  p <- p + scale_fill_manual(values = c("#BE8AD1", "#97B0E5")) + theme(legend.position = 'none')
  p
}

b.vs.w.plot(merged.with.race, mode = "actual")
b.vs.w.plot(merged.with.race, mode = "naive")
b.vs.w.plot(merged.with.race, mode = "predicted")

# actual tmb < 50
b.vs.w.plot(merged.with.race[actual.tmb  < 50], mode = "actual") + ylim(0, 60) + theme.figure
ggsave('b.vs.w.actual.pdf', width = 2, height = 3)
b.vs.w.plot(merged.with.race[actual.tmb < 50], mode = "naive") + ylim(0, 60) + theme.figure
ggsave('b.vs.w.naive.pdf', width = 2, height = 3)
b.vs.w.plot(merged.with.race[actual.tmb < 50], mode = "predicted") + ylim(0, 60) + theme.figure
ggsave('b.vs.w.predicted.pdf', width = 2, height = 3)


wilcox.test(merged.with.race[actual.tmb < 50][race.f == 'White']$naive.tmb, merged.with.race[actual.tmb < 50][race.f == 'Black']$naive.tmb)
#p-value = 3.939e-07
wilcox.test(merged.with.race[actual.tmb < 50][race.f == 'White']$actual.tmb, merged.with.race[actual.tmb < 50][race.f == 'Black']$actual.tmb)
# n.s.
wilcox.test(merged.with.race[actual.tmb < 50][race.f == 'White']$predicted.tmb, merged.with.race[actual.tmb < 50][race.f == 'Black']$predicted.tmb)
# p-value = 0.0001

median.predicted.w <- median(merged.with.race[actual.tmb < 50][race.f == 'White']$predicted.tmb)
# 1.85
median.predicted.b <- median(merged.with.race[actual.tmb < 50][race.f == 'Black']$predicted.tmb)
# 3.439
3.439/1.85
# 1.8 times higher.
median.naive.w <- median(merged.with.race[actual.tmb < 50][race.f == 'White']$naive.tmb)
# 11.146
median.naive.b <- median(merged.with.race[actual.tmb < 50][race.f == 'Black']$naive.tmb)
# 30.63
30.63/11.146 
# 2.75 times higher.
median.actual.w <- median(merged.with.race[actual.tmb < 50][race.f == 'White']$actual.tmb)
# 11.146
median.actual.b <- median(merged.with.race[actual.tmb < 50][race.f == 'Black']$actual.tmb)

# does PureCN eliminate bias??
full.venn.vars <- fread('../vs_purecn/tcga/pure_tabnet_merge.csv')
# get rid of train patients
test.val.venn.vars <- full.venn.vars[sample %in% data.splits[split %in% c('test','validation')]$patient]

merged.with.race <- merge(test.val.venn.vars, tcga.clin[,.(bcr_patient_barcode, race, ethnicity, jewish_origin)], by.x = 'sample', by.y = 'bcr_patient_barcode')

merged.with.race[, race.f := 'Unknown']
merged.with.race[race == "WHITE" , race.f := 'White']
merged.with.race[race == "BLACK OR AFRICAN AMERICAN" , race.f := 'Black']

true.tmb <- merged.with.race[truth == 1]
true.tmb[,true.tmb := .N, by = sample]

non.outlier.patients <- unique(true.tmb[true.tmb < (50*41)]$sample)

somatic.purecn.race <- merged.with.race[POSTERIOR.SOMATIC > 0.005][sample %in% non.outlier.patients]

merged.with.race[POSTERIOR.SOMATIC > 0.005] # 23 thousand rows.
merged.with.race[tabnet_proba_1 > 0.5] # 38 thousand rows.

tmb.purecn.race <- somatic.purecn.race[, purecn.tmb := .N, by = sample]

purecn.dt <- unique(tmb.purecn.race[,.(purecn.tmb, race.f, sample)])

ggplot(purecn.dt[race.f %in% c("Black", "White")], aes(x = race.f, fill = race.f))  + xlab('Race') + geom_boxplot(aes(y = purecn.tmb)) + ylab('PureCN Tumor-Only TMB') 

median(purecn.dt[race.f == 'White']$purecn.tmb)
# 50, so tmb = 1.22

median(purecn.dt[race.f == 'Black']$purecn.tmb)
# 99, so tmb = 2.41.  
2.41 / 1.22
# 1.98 fold difference