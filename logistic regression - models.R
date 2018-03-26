#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Using the primary dataset and the list of candidate metabolites as inputs, generate
#' regression models, predicted probabilities and actual outcomes for ROC curves.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/metabolomics.relapse.v20180306.1.rdata')

#' For reference: lists of the selected metabolites for each outcome.
#' mrd.metabs <- c('pyruvate', "`cytidine 5'-diphosphocholine`", '`1-arachidonylglycerol (20:4)`')
#' relapse.metabs <- c('`1-linoleoyl-GPI (18:2)*`', '`10-nonadecenoate (19:1n9)`', 'valine' , '`gamma-CEHC`')



# MRD models --------------------------------------------------------------

mrd.clin.model <- glm(mrd ~ nci.risk + immunophenotype + cyto.two.cat, data = subset(met, !is.na(met$mrd)), 
                      family = binomial(link = 'logit'))

mrd.metab.model <- glm(mrd ~ pyruvate + `cytidine 5'-diphosphocholine` + `1-arachidonylglycerol (20:4)`,
                       data = subset(met, !is.na(met$mrd)), family=binomial(link='logit'))

mrd.combined.model <- glm(mrd ~ pyruvate + `cytidine 5'-diphosphocholine` + `1-arachidonylglycerol (20:4)` + nci.risk + immunophenotype + cyto.two.cat,
                          data = subset(met, !is.na(met$mrd)), family=binomial(link='logit'))

#' Elements needed for ROC curve analysis.
mrd.outcomes <- ifelse(as.numeric(mrd.clin.model$data$mrd)==2,1,0) #' Vector of outcomes will be the same for all 3 models.
mrd.clin.model.predprob <- mrd.clin.model$fitted.values
mrd.metab.model.predprob <- mrd.metab.model$fitted.values
mrd.combined.model.predprob <- mrd.combined.model$fitted.values

mrd.models <- list(clin.model = list(model = mrd.clin.model, outcomes = mrd.outcomes, predicted.probabilities = mrd.clin.model.predprob),
                   metab.model = list(model = mrd.metab.model, outcomes = mrd.outcomes, predicted.probabilities = mrd.metab.model.predprob),
                   combined.model = list(model = mrd.combined.model, outcomes = mrd.outcomes, predicted.probabilities = mrd.combined.model.predprob))

#' Methods for accessing stored elements of this list:
#' tmp <- mrd.models$clin.model #' First nested list
#' tmp$clin.model$predicted.probabilities #' Specified element of first nested list

save(mrd.models, file = './Datasets/Expanded datasets/Logistic regression models/logreg.models.mrd.v20180326.1.rdata')



# Relapse models ----------------------------------------------------------

rel.clin.model <- glm(relapse ~ nci.risk + immunophenotype + cyto.two.cat + mrd, data = subset(met, !is.na(met$relapse) & !is.na(met$mrd)), 
                      family = binomial(link = 'logit'))

rel.metab.model <- glm(relapse ~ `1-linoleoyl-GPI (18:2)*` + `10-nonadecenoate (19:1n9)` + valine + `gamma-CEHC`,
                       data = subset(met, !is.na(met$relapse) & !is.na(met$mrd)), family=binomial(link='logit'))

rel.combined.model <- glm(relapse ~ `1-linoleoyl-GPI (18:2)*` + `10-nonadecenoate (19:1n9)` + valine + `gamma-CEHC` + nci.risk + immunophenotype + cyto.two.cat + mrd,
                          data = subset(met, !is.na(met$relapse) & !is.na(met$mrd)), family=binomial(link='logit'))

rel.outcomes <- ifelse(as.numeric(rel.clin.model$data$relapse)==2,1,0) #' Vector of outcomes will be the same for all 3 models.
rel.clin.model.predprob <- rel.clin.model$fitted.values
rel.metab.model.predprob <- rel.metab.model$fitted.values
rel.combined.model.predprob <- rel.combined.model$fitted.values

rel.models <- list(clin.model = list(model = rel.clin.model, outcomes = rel.outcomes, predicted.probabilities = rel.clin.model.predprob),
                   metab.model = list(model = rel.metab.model, outcomes = rel.outcomes, predicted.probabilities = rel.metab.model.predprob),
                   combined.model = list(model = rel.combined.model, outcomes = rel.outcomes, predicted.probabilities = rel.combined.model.predprob))

save(rel.models, file = './Datasets/Expanded datasets/Logistic regression models/logreg.models.relapse.v20180326.2.rdata')

rm(list = ls()); gc()
