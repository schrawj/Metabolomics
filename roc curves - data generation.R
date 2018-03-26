#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Generate ROC objects for model comparison and plotting.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')
load('./Expanded datasets/Logistic regression models/logreg.models.mrd.v20180326.1.rdata')
load('./Expanded datasets/Logistic regression models/logreg.models.relapse.v20180326.2.rdata')

require(pROC)



# Generate ROC objects: MRD -----------------------------------------------

mrd.roc.clin <- roc(response = mrd.models$clin.model$outcomes, predictor = mrd.models$clin.model$predicted.probabilities)
mrd.roc.metab <- roc(response = mrd.models$metab.model$outcomes, predictor = mrd.models$metab.model$predicted.probabilities)
mrd.roc.combined <- roc(response = mrd.models$combined.model$outcomes, predictor = mrd.models$combined.model$predicted.probabilities)

mrd.roc.objects <- list(clinical = mrd.roc.clin, metabolite = mrd.roc.metab, combined = mrd.roc.combined)

save(mrd.roc.objects, file = './Expanded datasets/ROC objects/mrd.roc.objects.rdata')



# Generate ROC objects: relapse -------------------------------------------

rel.roc.clin <- roc(response = rel.models$clin.model$outcomes, predictor = rel.models$clin.model$predicted.probabilities)
rel.roc.metab <- roc(response = rel.models$metab.model$outcomes, predictor = rel.models$metab.model$predicted.probabilities)
rel.roc.combined <- roc(response = rel.models$combined.model$outcomes, predictor = rel.models$combined.model$predicted.probabilities)

rel.roc.objects <- list(clinical = rel.roc.clin, metabolite = rel.roc.metab, combined = rel.roc.combined)

save(rel.roc.objects, file = './Expanded datasets/ROC objects/rel.roc.objects.rdata')

rm(list = ls()); gc()