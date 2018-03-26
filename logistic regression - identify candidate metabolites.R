#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Identify metabolites to include in logistic regression models for MRD and relapse.
#' 
#' Begin with all that cleared univariate significance and choose a subset by backwards
#' selection.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



# Load in metabolite significance data ------------------------------------

require(dplyr)

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.2.rdata')
load('./Datasets/metabolomics.relapse.v20180306.1.rdata')



# Choose metabolites: MRD -------------------------------------------------

mrd.can <- filter(met.signif, mrd.pvalue.kruskal <= 0.001 | mrd.pvalue.t.test <= 0.001)
mrd.can <- c('mrd', mrd.can$compound)

#' Subset main dataset.
#' Keep rows with non-missing MRD values and columns for candidate metabolites.
mrd <- met[!is.na(met$mrd) , names(met) %in% mrd.can]

#' Baseline model: all candidate metabolites.
mrd.model <- glm(mrd ~ ., data = mrd, family = binomial(link = 'logit'))

#' Stepwise backwards selection based on model AIC.
#' Top 5 retained metabolites are - in order of greatest to least AIC 
#' increase upon exclusion - pyruvate, cytidine 5'-disphosphocholine, 
#' 1-arachidonylglycerol (20:4), 2-hydroxyoctanoate and phsophate.
mrd.back <- step(mrd.model)

curated.mrd.model <- glm(mrd ~ pyruvate + `cytidine 5'-diphosphocholine` + `1-arachidonylglycerol (20:4)`,
                         data = mrd, family=binomial(link='logit'))



# Choose metabolites: relapse ---------------------------------------------

#' More relaxed criterion for relapse.
rel.can <- filter(met.signif, -log10(relapse.pvalue.kruskal) >= 2.5 | -log10(relapse.pvalue.t.test) >= 2.5)
rel.can <- c('relapse', rel.can$compound)

rel <- met[!is.na(met$relapse), names(met) %in% rel.can]

rel.model <- glm(relapse ~ ., data = rel, family = binomial(link = 'logit'))

#' Stepwise backwards selection based on model AIC.
#' Top 5 retained metabolites are - in order of greatest to least AIC 
#' increase upon exclusion - 1-linoleoyl-GPI (18:2*), 
#' 10-nonadecenoate (19:1n9), valine and gamma-CEHC.
rel.back <- step(rel.model)

curated.relapse.model <- glm(relapse ~ `1-linoleoyl-GPI (18:2)*` + `10-nonadecenoate (19:1n9)` + valine + `gamma-CEHC`,
                             data = rel, family=binomial(link='logit'))



# Generate and save a list of top metabolites -----------------------------

top.metabs <- list(mrd.metabs = c('pyruvate', "`cytidine 5'-diphosphocholine`", '`1-arachidonylglycerol (20:4)`'),
                   relapse.metabs = c('`1-linoleoyl-GPI (18:2)*`', '`10-nonadecenoate (19:1n9)`', 'valine' , '`gamma-CEHC`'))

save(top.metabs, file = './Datasets/Expanded datasets/metabolites.for.lr.models.20180326.1.rdata')

rm(list = ls()); gc()