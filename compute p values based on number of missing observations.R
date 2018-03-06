#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.06.14.
#' 
#' We know that for some of the less commonly detected compounds, t tests do not show 
#' any statistical significance despite large fold changes.  For an even smaller 
#' number of compounds, the significance of the findings can be increased by computing 
#' a dummy variable, which equals 1 if the compound was detected and 0 if it was not.
#' 
#' Split the compounds into two camps: thsoe detected in more than 50 samples, and those
#' detected in 50 or less.
#' 
#' For compounds detected in > 50 samples, compute p value by t test.
#' 
#' For compounds detected in =< samples, dummy code and compute p value by Fisher's 
#' exact test.
#' 
#' Bind all compounds back into a single dataset and re-do the Volcano plots.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' Gives us the number of NAs.
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.compounds.expanded.info.v20170614.rdata")

#' For the metabolites we want to treat continuously.
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20170612.2.rdata")

#' For the metabolites we want to treat categorically.
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolites.original.scale.v20170531.1.rdata")

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Cateogrical compounds.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' It's missing outcome data.
met.orig <- left_join(metab.orig,
                      
                      select(met, id, mrd, relapse),
                      
                      by= 'id') 

#' It problematically contains some drug metabolites that are missing in 99 samples.
met.orig <- select(met.orig, -`2-phosphoglycerate`, -`gabapentin`, -`hydroquinone beta-D-glucopyranoside`,
                   -`hydroxypioglitazone (M-IV)`, -`dehydrolithocholate`, -`ibuprofen acyl glucuronide`,
                   -`pantoprazole `, -`pioglitazone`, -`quetiapine`, -`X - 24585`)


met.orig <- met.orig[,c(1,969,970,2:968)]

save(met.orig, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.compounds.original.scale.v20170614.2.rdata')

#' Exlucde compounds missing in exactly 99 patients.
metabolite.annotations <- filter(metabolite.annotations, na.counts < 99)
#' A vector of columns we wish to keep.
missing50 <- subset(metabolite.annotations, na.counts >= 50)
missing50 <- c('id','mrd', 'relapse', missing50$compound)

#' A data frame containing measurements on the 178 compounds we will treat categorically.
test.categorical <- subset(met.orig[,missing50])

#' Print a few samples.
test.categorical[1:10,2:7]
test.categorical[51:60, 8:15]

#' A function to recode the exising variable.
dummyvar <- function(x){
                          ifelse(is.na(x),0,1)
                       }

#' Apply across the metabolite columns.
test.categorical[,4:181] <- apply(test.categorical[,4:181], 2, dummyvar)

#' Compute Fisher's exact test p values for dummy coded compounds.
categorical.compounds <- data.frame(compound = as.character(missing50[4:181]),
                                 
                                    mrd.pvalue = apply(test.categorical[,4:181], 2, function(x){
                                   
                                                  f <- with(test.categorical, fisher.test(x, mrd))
                                                  f <- f$p.value}),
                                 
                                    relapse.pvalue = apply(test.categorical[,4:181], 2, function(x){
                                   
                                          ff <- with(test.categorical, fisher.test(x, relapse))
                                          ff <- ff$p.value}))

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Continuous compounds.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' A vector of columns we wish to keep.
present50 <- subset(metabolite.annotations, na.counts < 50)
present50 <- c('id','mrd','relapse',present50$compound)

#' Should only contain id, MRD and relapse.
#' ...Which is exactly what it contains.
print(intersect(missing50, present50))

#' Remove unwanted columns.
test.continuous <- subset(met[,present50])

continuous.compounds <- data.frame(compound = as.character(present50[4:792]), 
                                  
                                  mrd.pvalue =
                                    
                                    apply(test.continuous[,c(4:792)], 2, function(x){
                                      
                                      t <- t.test(x ~ mrd, data = test.continuous, na.rm = TRUE)
                                      t <- t$p.value}),
                                  
                                  relapse.pvalue = 
                                    
                                    apply(test.continuous[,c(4:792)], 2, function(x){
                                      
                                      tt <- t.test(x ~ relapse, data = test.continuous, na.rm = TRUE)
                                      tt <- tt$p.value}))

#' Put them back together.
p.values.split.applied <- rbind(continuous.compounds, categorical.compounds)

#' Add the mean ratios we computed previously.

p.values.split.applied <- left_join(p.values.split.applied,
                                    
                                    select(met.multiob.pvalues, compound, ratio.of.means.mrd, ratio.of.means.relapse, na.counts),
                                    
                                    by = 'compound')

#' A flag for method used to obtain p value.
p.values.split.applied$method <- as.character(ifelse(p.values.split.applied$na.counts >= 50,
                                                     'p-value by Chi Squared test',
                                                     'p-value by welch two-sample t test'))

save(p.values.split.applied, 
     file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/p.values.and.mean.ratios.split.method.rdata')

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

