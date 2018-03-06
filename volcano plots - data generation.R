#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2018.06.03.
#' 
#' Generate univariate p-values for metabolites in connection with MRD and relapse in the 
#' newest dataset with N=99 subjects and more recent endpoints.
#' 
#' P-values will be computed using Student's t-test and Kruskal-Wallis test.
#' 
#' Many of the previous tasks in this script, e.g. filtering out compounds detected in 
#' just a single sample, have already been accomplished and have been removed.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



# Generate p values -------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/metabolomics.relapse.v20180306.1.rdata')

#' Generates p-values for each compound and each endpoint according to 
#' one parametric and one non.parametric method.
#' Also generates fold change data needed for subsequent volcano plots.
met.signif <- data.frame(compound = names(met[35:746]), 
                     
                     fold.change.mrd =                        
                       
                       apply(met[,35:746], 2, function(x){
                           tab <- aggregate(x ~ mrd, data = met, mean)
                           tab[2,2]/tab[1,2]}),
                                     
                     mrd.pvalue.t.test =
                                     
                       apply(met[,35:746], 2, function(x){
                           t <- t.test(x ~ mrd, data = met, na.rm = TRUE)
                           t <- t$p.value}),
                     
                     mrd.pvalue.kruskal = 
                       
                       apply(met[,35:746], 2, function(x){
                           k <- kruskal.test(x ~ mrd, data = met)
                           k <- k$p.value}),
                     
                     fold.change.relapse =                        
                       
                       apply(met[,35:746], 2, function(x){
                           tab <- aggregate(x ~ relapse, data = met, mean)
                           tab[2,2]/tab[1,2]}),
                       
                     relapse.pvalue.t.test = 
                                     
                         apply(met[,35:746], 2, function(x){
                           tt <- t.test(x ~ relapse, data = met, na.rm = TRUE)
                           tt <- tt$p.value}),
                     
                     relapse.pvalue.kruskal = 
                       
                         apply(met[,35:746], 2, function(x){
                           k <- kruskal.test(x ~ relapse, data = met)
                           k <- k$p.value}))

row.names(met.signif) <- 1:712

save(met.signif, file = './Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.1.rdata')

rm(met.signif, met); gc()



# Merge with compound metadata --------------------------------------------

require(dplyr)

#' Certain useful information on biological pathways and missingness was
#' previously compiled.  Join to new significance measures data frame.
setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/Expanded datasets/metabolomics.compounds.expanded.info.v20170803.rdata')
load('./Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.1.rdata')

met.signif <- left_join(met.signif,
                        select(metabolite.annotations, compound, super.pathway, sub.pathway, na.counts, KEGG, HMDB),
                        by = 'compound')

save(met.signif, file = './Datasets/Expanded datasets/metabolomics.compounds.significance.v20180306.2.rdata')

rm(met.signif, metabolite.annotations); gc()

