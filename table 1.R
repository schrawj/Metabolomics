#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2018.03.06.
#' 
#' Table 1 for Metabolomics-Relapse manuscript.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------




# Prep environment --------------------------------------------------------

require(dplyr); require(gmodels)

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/metabolomics.relapse.v20180306.1.rdata')



# Generate results --------------------------------------------------------

table.vars <- c('sex','race.eth.collapsed','overwt.status','nci.risk','immunophenotype','cyto.two.cat')

for (i in table.vars){
  print(i)
  CrossTable(met[,i], met$mrd, prop.chisq = FALSE, prop.r = FALSE, prop.t = FALSE, chisq = TRUE, fisher = TRUE)
  CrossTable(met[,i], met$relapse, prop.chisq = FALSE, prop.r = FALSE, prop.t = FALSE, chisq = TRUE, fisher = TRUE)
}

endpoints <- c('mrd','relapse')

for (i in endpoints){
  print(aggregate(age.dx ~ met[,i], data = met, mean))
  print(aggregate(age.dx ~ met[,i], data = met, sd))
  print(kruskal.test(age.dx ~ met[,i], data = met))
}