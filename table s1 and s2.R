
# Generate tables S1 and S2 for Blood manuscript --------------------------

load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.2.rdata")

require(dplyr)

tmp <- rbind(filter(met.signif, mrd.pvalue.kruskal <= 0.001),
             filter(met.signif, mrd.pvalue.t.test <= 0.001))
tmp <- tmp[!duplicated(tmp$compound), ]

write.csv(select(tmp, compound, fold.change.mrd, mrd.pvalue.kruskal, mrd.pvalue.t.test),
          file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Tables and figures/Table S1.csv',
          row.names = FALSE)

tmp <- rbind(filter(met.signif, relapse.pvalue.kruskal <= 0.004),
             filter(met.signif, relapse.pvalue.t.test <= 0.004))
tmp <- tmp[!duplicated(tmp$compound), ]

write.csv(select(tmp, compound, fold.change.relapse, relapse.pvalue.kruskal, relapse.pvalue.t.test),
          file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Tables and figures/Table S2.csv',
          row.names = FALSE) 