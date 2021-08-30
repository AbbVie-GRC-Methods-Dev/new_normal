library(data.table)
library(ggplot2)

tcga.dir <- '../../test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/'
full.c.maf <- data.table::fread(paste0(tcga.dir, 'mafs/tabnet_trained_tcga_v3_predictions_all_features.csv'))

full.c.maf[,result := '']
full.c.maf[(truth == 1) & (tabnet_pred == 1), result := 'TP']
full.c.maf[(truth == 0) & (tabnet_pred == 1), result := 'FP']
full.c.maf[(truth == 0) & (tabnet_pred == 0), result := 'TN']
full.c.maf[(truth == 1) & (tabnet_pred == 0), result := 'FN']

perf.all <- data.table::fread('data_splits_210_with_tmb38.csv')

perf.all.no.laml <- perf.all[indication != 'LAML']

high.sens <- perf.all.no.laml[sens > median(perf.all.no.laml$sens)]$patient
low.sens <- perf.all.no.laml[sens < median(perf.all.no.laml$sens)]$patient

h.sens <- full.c.maf[sample %in% high.sens]
l.sens <- full.c.maf[sample %in% low.sens]

mean(h.sens[truth == 1]$t_alt_freq)
mean(l.sens[truth == 1]$t_alt_freq)

mean(full.c.maf[result == 'TP']$t_alt_freq)
mean(full.c.maf[result == 'FP']$t_alt_freq)
mean(full.c.maf[result == 'FN']$t_alt_freq)
