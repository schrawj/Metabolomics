#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Prepare input for IMPaLA enrichment analysis.
#' 
#' http://impala.molgen.mpg.de/
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------



# Prep environment --------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')
load("./metabolomics.relapse.v20180306.1.rdata")
load("./Expanded datasets/metabolomics.compounds.expanded.info.v20170803.rdata")



# Generate log2-scale mean concentration data by MRD status ---------------

require(dplyr)

impala <- data.frame(mrdpos = 1, mrdneg = 1)

for (i in met[ , 35:746]){
  tmp2 <- aggregate(i ~ mrd, data = met, mean)
  tmp2 <- data.frame(mrdpos = tmp2[2,2], mrdneg = tmp2[1,2])
  impala <- rbind(impala, tmp2)
}

impala <- impala[-1,]
impala <- transmute(impala, logmrdpos = log2(mrdpos), logmrdneg = log2(mrdneg))
impala$compound <- colnames(met[35:746])
impala <- left_join(impala, 
                 select(metabolite.annotations, compound, KEGG, KEGG.generic, HMDB), 
                 by = 'compound')
impala <- select(impala, -compound)



# Write table using KEGG ID -----------------------------------------------

require(dplyr)

#' KEGG or HMDB ID will produce generally equivalent results.
impala <- impala[grepl('^C', impala$KEGG.generic), ]
impala <- select(impala, KEGG.generic, logmrdpos, logmrdneg)

write.table(impala, file='C:/Users/schraw/Desktop/IMPaLA_KEGG.txt',row.names = FALSE, sep = '\t',
            quote = FALSE)



# Write table using HMDB ID -----------------------------------------------

require(dplyr)

impala <- impala[grepl('^H', impala$HMDB), ]
impala <- select(impala, HMDB, logmrdpos, logmrdneg)

write.table(impala, file='C:/Users/schraw/Desktop/IMPaLA_HMDB.txt',row.names = FALSE, sep = '\t',
            quote = FALSE)



# Generate log2-scale mean concentration data by relapse status -----------

require(dplyr)

impala <- data.frame(relapse.n = 1, relapse.y = 1)

for (i in met[ , 35:746]){
  tmp2 <- aggregate(i ~ relapse, data = met, mean)
  tmp2 <- data.frame(relapse.n = tmp2[1,2], relapse.y = tmp2[2,2])
  impala <- rbind(impala, tmp2)
}

impala <- impala[-1,]
impala <- transmute(impala, logrelapse.n = log2(relapse.n), logrelapse.y = log2(relapse.y))
impala$compound <- colnames(met[35:746])
impala <- left_join(impala, 
                    select(metabolite.annotations, compound, KEGG, KEGG.generic, HMDB), 
                    by = 'compound')
impala <- select(impala, -compound)



# Write table using KEGG ID -----------------------------------------------

require(dplyr)

#' KEGG or HMDB ID will produce generally equivalent results.
impala <- impala[grepl('^C', impala$KEGG.generic), ]
impala <- select(impala, KEGG.generic, logrelapse.n, logrelapse.y)

write.table(impala, file='C:/Users/schraw/Desktop/IMPaLA_KEGG.txt',row.names = FALSE, sep = '\t',
            quote = FALSE)



# Write table using HMDB ID -----------------------------------------------

require(dplyr)

impala <- impala[grepl('^H', impala$HMDB), ]
impala <- select(impala, HMDB, logrelapse.n, logrelapse.y)

write.table(impala, file='C:/Users/schraw/Desktop/IMPaLA_HMDB.txt',row.names = FALSE, sep = '\t',
            quote = FALSE)




# Sanity check ------------------------------------------------------------

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' A negative control of sorts.  Generate a data frame with the same compound IDs, but 
#' random normal variables.  Confirm that the same pathways are not enriched. 
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
ids <- metabolite.annotations[grepl('^H', metabolite.annotations$HMDB), ]
ids <- data.frame(id = ids$HMDB)

values <- data.frame(group1 = rnorm(604, 0, 1),
                     group2 = rnorm(604, 0, 1))

spurious <- cbind(ids, values)

write.table(spurious, file='C:/Users/schraw/Desktop/IMPaLA_spurious.txt',row.names = FALSE, sep = '\t',
            quote = FALSE)

