
# Generate tables S1 and S2 for Blood manuscript --------------------------

#' New approach to explaining candidate metabolites: filter currency metabolites, then keep top 20 by K-W p-value.
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.2.rdata")

require(dplyr)

tmp <- arrange(met.signif, mrd.pvalue.kruskal)
#' Relevant currency metabolites to exclude.
tmp <- filter(tmp, !(compound %in% c('phosphate', "cytidine 5'-diphosphocholine","inosine 5'-monophosphate (IMP)", "adenosine 3',5'-cyclic monophosphate (cAMP)", "adenosine 3'-monophosphate (3'-AMP)")))
tmp <- tmp[1:20, ]

write.csv(select(tmp, compound, fold.change.mrd, mrd.pvalue.kruskal),
          file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Tables and figures/Table S1 v20180426.csv',
          row.names = FALSE)

tmp <- arrange(met.signif, relapse.pvalue.kruskal)
#' Relevant currency metabolites to exclude.
tmp <- filter(tmp, !(compound %in% c('phospate', "inosine 5'-monophosphate (IMP)", "adenosine 3',5'-cyclic monophosphate (cAMP)", "adenosine 3'-monophosphate (3'-AMP)")))
tmp <- tmp[1:20, ]

write.csv(select(tmp, compound, fold.change.relapse, relapse.pvalue.kruskal),
          file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Tables and figures/Table S2 v20180426.csv',
          row.names = FALSE)