#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Evaluate and plot ROC objects.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')
load('./Expanded datasets/ROC objects/mrd.roc.objects.rdata')
load('./Expanded datasets/ROC objects/rel.roc.objects.rdata')

require(pROC)



# DeLong's test for correlated ROC curves ---------------------------------

roc.test(mrd.roc.objects$clinical, mrd.roc.objects$metabolite) # p=0.08595
roc.test(mrd.roc.objects$clinical, mrd.roc.objects$combined) # p=0.002202

roc.test(rel.roc.objects$clinical, rel.roc.objects$metabolite) # p=0.4586
roc.test(rel.roc.objects$clinical, rel.roc.objects$combined) # p=0.03966



# Plots -------------------------------------------------------------------

plot.roc(mrd.roc.objects$clinical, 
         print.auc = TRUE, max.auc.polygon = TRUE, print.auc.y = 0.35, print.auc.cex = 1.2)
plot.roc(mrd.roc.objects$metabolite, 
         print.auc = TRUE, print.auc.cex = 1.2, print.auc.y = 0.4, col = 'blue', add = TRUE)
plot.roc(mrd.roc.objects$combined, 
         print.auc = TRUE, print.auc.cex = 1.2, col = 'red', print.auc.y = 0.45, add = TRUE)

#' Caption can be added as below.  Leave blank for now; looks better if added in post-production.

#' legend('bottomright', c('Clinical variables', 'Metabolites, p = 0.09 vs. clinical model', 
#'                        'Combination, p = 0.002 vs. clinical model'), 
#'       text.col = c('black','blue','red'), cex = 1.2)

plot.roc(rel.roc.objects$clinical, 
         print.auc = TRUE, max.auc.polygon = TRUE, print.auc.y = 0.35, print.auc.cex = 1.2)
plot.roc(rel.roc.objects$metabolite, 
         print.auc = TRUE, print.auc.cex = 1.2, print.auc.y = 0.4, col = 'blue', add = TRUE)
plot.roc(rel.roc.objects$combined, 
         print.auc = TRUE, print.auc.cex = 1.2, col = 'red', print.auc.y = 0.45, add = TRUE)

#' Caption can be added as below.  Leave blank for now; looks better if added in post-production.

#' legend('bottomright', c('Clinical variables', 'Metabolites, p = 0.46 vs. clinical model', 
#'                        'Combination, p = 0.04 vs. clinical model'), 
#'       text.col = c('black','blue','red'), cex = 1.2)

rm(mrd.roc.objects, rel.roc.objects); gc()