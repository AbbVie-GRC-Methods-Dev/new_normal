
hugo.dir <- '../../test/no_ps_fit_with_snp-vaf-bin/hugo_melanoma/'
data.splits.h <- data.table::fread(paste0(hugo.dir,'data_splits.csv'))

full.h.maf  <- data.table::fread(paste0(hugo.dir, 'mafs/tabnet_trained_tcga_v3_predictions_all_features.csv'))


full.h.maf[,result := '']
full.h.maf[(truth == 1) & (tabnet_pred == 1), result := 'TP']
full.h.maf[(truth == 0) & (tabnet_pred == 1), result := 'FP']
full.h.maf[(truth == 0) & (tabnet_pred == 0), result := 'TN']
full.h.maf[(truth == 1) & (tabnet_pred == 0), result := 'FN']


data.splits.h[, sens :=2] 
data.splits.h[, spec :=2]
data.splits.h[, ppv :=2]
data.splits.h[, npv :=2]
data.splits.h[, f1 :=2]
# without normal, with tabnet.
data.splits.h[, predicted.tmb :=-1]
# with normal.
data.splits.h[, actual.tmb :=-1]

perf.h <- process.maf.subset(full.h.maf, data.splits.h)

perf.h.indels <- process.maf.subset(full.h.maf[nchar(ref) != nchar(alt)], data.splits.h)
perf.h.snvs <- process.maf.subset(full.h.maf[nchar(ref) == nchar(alt)], data.splits.h)

rename.metrics.1 <- function(input.dt){
   dt <- copy(input.dt)
   setnames(old = 'sens', new = 'TPR',  dt)
   setnames(old = 'spec', new = 'TNR',  dt)
   setnames(old = 'ppv', new = 'PPV',  dt)
   setnames(old = 'npv', new = 'NPV',  dt)
   return(dt)
}

renamed.perf.h <- rename.metrics.1(perf.h)

perf.h.m <- melt(renamed.perf.h, measure.vars = c('TPR','TNR','PPV','NPV'), variable.name = 'metric' )

ggplot(perf.h.m) + geom_boxplot(aes(x = metric, y = value, fill = metric )) + tcga_boxplot_theme + scale_fill_manual(values = four.beach.colors) + facet_grid(~split) + ylim(0, 1)
ggsave('hugo.metrics.pdf', height = 4, width = 2.5)

