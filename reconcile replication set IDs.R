#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.23.
#' 
#' I have a few new spreadsheets to review, which contain MRNs, names, MRD
#' and relapse info for sets of children who might have sample available
#' for our replication set.
#' 
#' 'additional patients 6-17-18 with vials available' has info for patients
#' ID'd by Karen, for whom Theresa Ty at RTSS has already verified sample 
#' availabillity.
#' 
#' 'MRD+ Relapse' has info for patients Olga and Sydney identified.  This 
#' is a long list, but it needs to be cross-checked against our original
#' patient set to make sure there are no duplicates.  Karen also suggested
#' that I need to verify that all children were originally DX'd with ALL
#' (there may be some AML or lymphoma cases).  After we have verified that
#' all these patients are eligible cases not in our original study sample,
#' Theresa Ty can verify sample availability.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

#'TODO: Remove children originally DX'd with any condition than ALL.
#'TODO: Return curated list to Karen, Olga, Theresa Ty.




# Read in patient lists, remove duplicates, write joined list -------------

require(xlsx); require(dplyr)

karen.list <- read.xlsx(file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/additional patients 6-17-18 with vials available.xlsx',
                        sheetName = "Rabin's result_1_3_18", colIndex = 1:9, rowIndex = 1:58, header = TRUE, stringsAsFactors = FALSE)

olga.list.mrd <- read.xlsx(file = "Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/MRD+Relapse.xlsx", header = FALSE, stringsAsFactors = FALSE,
                           sheetName = 'Sheet1', colIndex = 1:10, rowIndex = 2:49)

olga.list.relapse <- read.xlsx(file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/MRD+Relapse.xlsx', header = FALSE, stringsAsFactors = FALSE,
                               sheetIndex = 1, colIndex = 1:10, rowIndex = 53:92)

load('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.subjects.expanded.identifiers.rdata')

olga.names <- c('unknown','mrn','rtss','unknown2','last.name','first.name','unknown3','unknown.date','mrd','relapse')

names(metab.clin.expanded.ids) <- tolower(names(metab.clin.expanded.ids))
names(olga.list.mrd) <- olga.names
names(olga.list.relapse) <- olga.names
names(karen.list) <- tolower(names(karen.list))

#' Use anti-joins to filter the kids in our current dataset.
#' Use semi-joins to generate lists of duplicated IDs.
mrd.dups <- semi_join(olga.list.mrd, metab.clin.expanded.ids, 'mrn')
olga.list.mrd <- anti_join(olga.list.mrd, metab.clin.expanded.ids, by = 'mrn') 

relapse.dups <- semi_join(olga.list.relapse, metab.clin.expanded.ids, 'mrn')
olga.list.relapse <- anti_join(olga.list.relapse, metab.clin.expanded.ids, 'mrn')

karen.dups <- semi_join(karen.list, metab.clin.expanded.ids, 'mrn')
karen.list <- anti_join(karen.list, metab.clin.expanded.ids, 'mrn')

olga.dups <- semi_join(olga.list.mrd, olga.list.relapse, 'mrn')
olga.list.mrd <- anti_join(olga.list.mrd, olga.list.relapse, 'mrn') 
olga.list <- arrange(rbind(olga.list.mrd, olga.list.relapse), last.name)

olga.dups <- semi_join(olga.list, karen.list, 'mrn')
olga.list <- anti_join(olga.list, karen.list, 'mrn')

olga.list$source.list <- 'Olga'
karen.list$source.list <- 'Karen'

joined.list <- rbind(select(karen.list, mrn, last.name, first.name, source.list),
                     select(olga.list, mrn, last.name, first.name, source.list))

write.xlsx(joined.list, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/olga.and.karen.lists.join.xlsx', row.names = FALSE)

rm(list = ls()); gc()
