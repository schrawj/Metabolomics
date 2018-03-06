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



# Random old code destined for deprecation --------------------------------

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Plots.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/p.values.and.mean.ratios.split.method.rdata")
require(ggplot2)
require(ggrepel)

#' Volcano plot, MRD.
volc.mrd <- ggplot(data=p.values.split.applied, aes(x=log2(ratio.of.means.mrd), y = -log10(mrd.pvalue))) + geom_point() + 
  geom_text_repel(aes(label=ifelse(-log10(mrd.pvalue) > 3, compound, '')), size = 3) +
  ggtitle('Volcano Plot for MRD', subtitle = 'Mean ratio for MRD+ : MRD- patients on x-axis; Statistical significance by various methods on y-axis.')
print(volc.mrd)

#' Volcano plot, relapse.
volc.relapse <- ggplot(data=p.values.split.applied, aes(x=log2(ratio.of.means.relapse), y = -log10(relapse.pvalue))) + 
  geom_point() + 
  geom_text_repel(aes(label=ifelse(-log10(relapse.pvalue) > 2.5, compound, '')), size = 3) +
  ggtitle('Volcano Plot for Relapse', subtitle = 'Mean ratio for relapsed:non-relapsed patients on x-axis; Statistical significance by various methods on y-axis.')
print(volc.relapse)

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.06.16.
#' 
#' Might be nice to know what sort of pathways these are in.  Append the Metabolon-
#' provided pathway data.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

require(xlsx)
require(dplyr)
setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')

#' Unfortunately, I don't already have this info in any R object.
metabolite.pathways <- read.xlsx('./Metabolon raw data/Original Metabolon dataset.xlsx',
                                 sheetName = 'OrigScale', rowIndex = 13:990, colIndex = 2:4,
                                 header = TRUE, stringsAsFactors = FALSE)

names(metabolite.pathways) <- tolower(names(metabolite.pathways))
metabolite.pathways <- rename(metabolite.pathways, compound = biochemical)

#' There are 8 super pathways and 80 subpathways.
metabolite.annotations <- left_join(metabolite.annotations, metabolite.pathways, by = 'compound')

#' Re-organize.
metabolite.annotations <- metabolite.annotations[,c(10:12,1:9)]

save(metabolite.annotations, 
     file = './metabolomics.compounds.expanded.info.v20170616.rdata')

#' Append to p values.
p.values.split.applied <- left_join(p.values.split.applied,
                                    
                                    select(metabolite.annotations, compound, super.pathway, sub.pathway),
                                    
                                    by = 'compound')

save(p.values.split.applied, 
     file = './Expanded datasets/p.values.and.mean.ratios.split.method.v20170616.rdata')

#' Volcano plots with color indicating super pathway.
volc.mrd <- ggplot(data=p.values.split.applied, aes(x=log2(ratio.of.means.mrd), y = -log10(mrd.pvalue), color = super.pathway)) + 
  geom_point() + 
  geom_text_repel(aes(label=ifelse(-log10(mrd.pvalue) > 3, compound, '')), size = 3) +
  xlab(label = 'Log2 ratio of mean concentration in MRD+:MRD- plasma') +
  ylab(label = '-Log10 p-value') +
  ggtitle('Candidate Metabolites for MRD Models', subtitle = 'Mean ratio for MRD+ : MRD- patients on x-axis; Statistical significance by various methods on y-axis.')

volc.relapse <- ggplot(data=p.values.split.applied, aes(x=log2(ratio.of.means.relapse), y = -log10(relapse.pvalue), color = super.pathway)) + 
  geom_point() + 
  geom_text_repel(aes(label=ifelse(-log10(relapse.pvalue) > 2.4, compound, '')), size = 4) +
  xlab(label = 'Log2 ratio of mean concentration in relapsed:non-relapsed plasma') +
  ylab(label = '-Log10 p-value') +
  ggtitle('Candidate Metabolites for Relapse Models')

print(volc.mrd)
print(volc.relapse)

