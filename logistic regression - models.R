#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Using the primary dataset and the list of candidate metabolites as inputs, generate
#' regression models, predicted probabilities and actual outcomes for ROC curves.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/metabolomics.relapse.v20180306.1.rdata')



# MRD models --------------------------------------------------------------

mrd.clin.model <- glm(mrd ~ nci.risk + immunophenotype + cyto.two.cat, data = subset(met, !is.na(met$mrd)), 
                      family = binomial(link = 'logit'))

mrd.metab.model <- glm(mrd ~ pyruvate + `2-hydroxyoctanoate` + `1-arachidonylglycerol (20:4)`,
                       data = subset(met, !is.na(met$mrd)), family=binomial(link='logit'))

#' Originally chose cytidine 5'-diphosphocholine as one of the metabolites. 
#' Removed it because it was causing problems with linear separation,
#' most likely to due to large number of samples with imputed values.
#' Replacing with 2-hydroxyoctanoate solves problem, and produces better AIC.
mrd.combined.model <- glm(mrd ~ pyruvate + `2-hydroxyoctanoate` + `1-arachidonylglycerol (20:4)` + nci.risk + immunophenotype + cyto.two.cat,
                          data = subset(met, !is.na(met$mrd)), family=binomial(link='logit'))
mrd.backselect.model <- step(mrd.combined.model)

#' Elements needed for ROC curve analysis.
mrd.outcomes <- ifelse(as.numeric(mrd.clin.model$data$mrd)==2,1,0) #' Vector of outcomes will be the same for all 3 models.
mrd.clin.model.predprob <- mrd.clin.model$fitted.values
mrd.metab.model.predprob <- mrd.metab.model$fitted.values
mrd.combined.model.predprob <- mrd.combined.model$fitted.values

mrd.models <- list(clin.model = list(model = mrd.clin.model, outcomes = mrd.outcomes, 
                                     predicted.probabilities = mrd.clin.model.predprob, 
                                     backxform.coef = exp(cbind(mrd.clin.model$coefficients, confint(mrd.clin.model, level = 0.95)))),
                   metab.model = list(model = mrd.metab.model, outcomes = mrd.outcomes, 
                                      predicted.probabilities = mrd.metab.model.predprob,
                                      backxform.coef = exp(cbind(mrd.metab.model$coefficients, confint(mrd.metab.model, level = 0.95)))),
                   combined.model = list(model = mrd.combined.model, outcomes = mrd.outcomes, 
                                         predicted.probabilities = mrd.combined.model.predprob,
                                         backxform.coef = exp(cbind(mrd.combined.model$coefficients, confint(mrd.combined.model, level = 0.95)))),
                   backseletion.model = list(model = mrd.backselect.model, outcomes = mrd.outcomes, 
                                             predicted.probabilities = mrd.backselect.model$fitted.values,
                                             backxform.coef = exp(cbind(mrd.backselect.model$coefficients, confint(mrd.backselect.model, level = 0.95)))))

#' Methods for accessing stored elements of this list:
#' tmp <- mrd.models$clin.model #' First nested list
#' tmp$clin.model$predicted.probabilities #' Specified element of first nested list

save(mrd.models, file = './Datasets/Expanded datasets/Logistic regression models/logreg.models.mrd.v20180328.1.rdata')



# Relapse models ----------------------------------------------------------

rel.clin.model <- glm(relapse ~ nci.risk + immunophenotype + cyto.two.cat + mrd, data = subset(met, !is.na(met$relapse) & !is.na(met$mrd)), 
                      family = binomial(link = 'logit'))

rel.metab.model <- glm(relapse ~ `1-linoleoyl-GPI (18:2)*` + `10-nonadecenoate (19:1n9)` + valine + picolinate,
                       data = subset(met, !is.na(met$relapse) & !is.na(met$mrd)), family=binomial(link='logit'))

rel.combined.model <- glm(relapse ~ `1-linoleoyl-GPI (18:2)*` + `10-nonadecenoate (19:1n9)` + valine + picolinate + nci.risk + immunophenotype + cyto.two.cat + mrd,
                          data = subset(met, !is.na(met$relapse) & !is.na(met$mrd)), family=binomial(link='logit'))
rel.backselect.model <- step(rel.combined.model)

rel.outcomes <- ifelse(as.numeric(rel.clin.model$data$relapse)==2,1,0) #' Vector of outcomes will be the same for all 3 models.
rel.clin.model.predprob <- rel.clin.model$fitted.values
rel.metab.model.predprob <- rel.metab.model$fitted.values
rel.combined.model.predprob <- rel.combined.model$fitted.values

rel.models <- list(clin.model = list(model = rel.clin.model, outcomes = rel.outcomes, predicted.probabilities = rel.clin.model.predprob,
                                     backxform.coef = exp(cbind(rel.clin.model$coefficients, confint(rel.clin.model, level = 0.95)))),
                   metab.model = list(model = rel.metab.model, outcomes = rel.outcomes, predicted.probabilities = rel.metab.model.predprob,
                                      backxform.coef = exp(cbind(rel.metab.model$coefficients, confint(rel.metab.model, level = 0.95)))),
                   combined.model = list(model = rel.combined.model, outcomes = rel.outcomes, predicted.probabilities = rel.combined.model.predprob,
                                         backxform.coef = exp(cbind(rel.combined.model$coefficients, confint(rel.combined.model, level = 0.95)))),
                   backseletion.model = list(model = rel.backselect.model, outcomes = rel.outcomes, 
                                             predicted.probabilities = rel.backselect.model$fitted.values,
                                             backxform.coef = exp(cbind(rel.backselect.model$coefficients, confint(rel.backselect.model, level = 0.95)))))

#' Model outputs feed into the script 'roc curves - data generation'.
save(rel.models, file = './Datasets/Expanded datasets/Logistic regression models/logreg.models.relapse.v20180426.1.rdata')

rm(list = ls()); gc()
