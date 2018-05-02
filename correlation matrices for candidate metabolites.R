#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.04.27.
#' 
#' Generate and export correlation matrices for the top 20 compounds 
#' associated with MRD and relapse.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(dplyr)

load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20180306.1.rdata")
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.2.rdata")

tmp <- arrange(met.signif, mrd.pvalue.kruskal)
tmp <- filter(tmp, !(compound %in% c('phosphate', "cytidine 5'-diphosphocholine","inosine 5'-monophosphate (IMP)", "adenosine 3',5'-cyclic monophosphate (cAMP)", "adenosine 3'-monophosphate (3'-AMP)")))
tmp <- tmp[1:20, ]
mrd.compounds <- c(tmp$compound)

tmp <- arrange(met.signif, relapse.pvalue.kruskal)
#' Relevant currency metabolites to exclude.
tmp <- filter(tmp, !(compound %in% c('phospate', "inosine 5'-monophosphate (IMP)", "adenosine 3',5'-cyclic monophosphate (cAMP)", "adenosine 3'-monophosphate (3'-AMP)")))
tmp <- tmp[1:20, ]
rel.compounds <- c(tmp$compound)

rm(met.signif, tmp)

mrd.cor <- met[ , colnames(met) %in% mrd.compounds]
rel.cor <- met[ , colnames(met) %in% rel.compounds]

tmp <- cor(mrd.cor, method = 'spearman')
write.csv(tmp, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/R outputs/cor.matrix.mrd.compounds.csv')

tmp <- cor(rel.cor, method = 'spearman')
write.csv(tmp, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/R outputs/cor.matrix.relapse.compounds.csv')

rm(list = ls()); gc()