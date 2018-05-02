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

#' Candidate metabolites have already been identified.  See script "tables S1 and S2" for details.
mrd.compounds <- read.csv(file = './Tables and figures/Table S1 v20180426.csv', stringsAsFactors = FALSE)
mrd.compounds <- c('mrd', mrd.compounds$compound)

rel.compounds <- read.csv(file = './Tables and figures/Table S2 v20180426.csv', stringsAsFactors = FALSE)
rel.compounds <- c('relapse', rel.compounds$compound)



# Choose metabolites: MRD -------------------------------------------------

#' Subset main dataset.
#' Keep rows with non-missing MRD values and columns for candidate metabolites.
mrd <- met[!is.na(met$mrd) , names(met) %in% mrd.compounds]

#' Baseline model: all candidate metabolites.
mrd.model <- glm(mrd ~ ., data = mrd, family = binomial(link = 'logit'))

#' Stepwise backwards selection based on model AIC.
mrd.back <- step(mrd.model)



# Choose metabolites: relapse ---------------------------------------------

rel <- met[!is.na(met$relapse), names(met) %in% rel.compounds]

rel.model <- glm(relapse ~ ., data = rel, family = binomial(link = 'logit'))

#' Stepwise backwards selection based on model AIC.
#' Top 4 retained metabolites are - in order of greatest to least AIC 
#' increase upon exclusion - 1-linoleoyl-GPI (18:2*), 
#' 10-nonadecenoate (19:1n9), valine and picolinate.
rel.back <- step(rel.model)



# Generate and save a list of top metabolites -----------------------------

top.metabs <- list(mrd.metabs = c('pyruvate', "`cytidine 5'-diphosphocholine`", '`1-arachidonylglycerol (20:4)`'),
                   relapse.metabs = c('`1-linoleoyl-GPI (18:2)*`', '`10-nonadecenoate (19:1n9)`', 'valine' , 'picolinate'))

#' These metabolites feed into the script 'logistic regression - models'.
save(top.metabs, file = './Datasets/Expanded datasets/metabolites.for.lr.models.20180426.1.rdata')

rm(list = ls()); gc()