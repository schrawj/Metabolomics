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



# Read in patient lists, remove duplicates, write joined list -------------

require(xlsx); require(dplyr); require(stringr)

karen.list <- read.xlsx(file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/additional patients 6-17-18 with vials available.xlsx',
                        sheetName = "Rabin's result_1_3_18", colIndex = 1:9, rowIndex = 1:58, header = TRUE, stringsAsFactors = FALSE)

olga.list.mrd <- read.xlsx(file = "Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/MRD+Relapse.xlsx", header = FALSE, stringsAsFactors = FALSE,
                           sheetName = 'Sheet1', colIndex = 1:10, rowIndex = 2:49)

olga.list.relapse <- read.xlsx(file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/MRD+Relapse.xlsx', header = FALSE, stringsAsFactors = FALSE,
                               sheetIndex = 1, colIndex = 1:10, rowIndex = 53:92)

load('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.subjects.expanded.identifiers.rdata')

olga.names <- c('unknown','mrn','rtss','unknown2','last.name','first.name','unknown3','dob','mrd','relapse')

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

#' Export Olga's list to Excel for some manual annotation.
#' Will fill in a more standardized variable for MRD, and flag cases who might not be suitable:
#' MRD not measured on day 29, relapsed after BMT, transferred care, AML or Burkitt lymphoma
write.xlsx(olga.list, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/MRD+ Relapse_jms.xlsx', row.names = FALSE)

#' Read the edited version back in, remove unnecessary columns and clean up relapse data.
olga.list <- read.xlsx(file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/MRD+ Relapse_jms_edited.xlsx',
                       stringsAsFactors = FALSE, sheetIndex = 1)
olga.list <- select(filter(olga.list, exclude.flag != 1), -exclude.flag, -unknown3)
olga.list$mrd.status <- paste('mrd',olga.list$mrd.status)
olga.list$relapse <- str_trim(olga.list$relapse, side = 'both')
olga.list$relapse <- tolower(olga.list$relapse)
olga.list$relapse.status <- ifelse(olga.list$relapse == 'no relapse', "no relapse", "relapse")
olga.list$group <- paste(paste(olga.list$mrd.status, olga.list$relapse.status))

combinations <- unique(paste(olga.list$mrd.status, olga.list$relapse.status))

for (i in 1:length(combinations)){
  
  name <- combinations[i]
  
  tmp <- select(filter(olga.list, group == combinations[i]), -group)
  
  write.xlsx(tmp, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/replication.set.potential.ids.xlsx',
             row.names = FALSE, append = TRUE, sheetName = name)
}

rm(list = ls()); gc()



# Further checks to Olga's list -------------------------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.07.27.
#' 
#' Per email exchange with Karen on 7/25, some of the patients in this list
#' are AML.  Need to further check them to confirm their original DX was 
#' ALL, make sure that % MRD is for day 29, and that we retrieve date the 
#' date of their original diagnosis.
#' 
#' As of today, have combed through EPIC to get this data and there were 
#' indeed 11 AML patients and 1 biphenotypic (T/myeloid) patient.  Also a 
#' few patients who did not undergo all of their therapy at TXCH.
#' 
#' Made a command decision to exclude 1 who was lost to follow-up.  Will 
#' hold off on excluding patients who did not have their initial therapy
#' at TXCH even though their D29 marrows will not be available. 
#' 
#' Import verified data and sort by MRD and relapse status.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(xlsx); require(dplyr)

#' Read in verified data.
curated.list <- read.xlsx('Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/replication.set.potential.ids.epic.check.xlsx',
                          header = TRUE, stringsAsFactors = FALSE, sheetName = 'received')

#' Exclude AML and biphenotypic patients - also 1 patient lost to follow-up.
curated.list <- select(filter(curated.list, exclude.flag != 1), -exclude.flag)

write.xlsx(curated.list, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/replication.set.potential.ids.epic.check.xlsx',
           row.names = FALSE, append = TRUE, sheetName = 'final', showNA = FALSE)

rm(list = ls()); gc()



# Check additional new patients against existing --------------------------

#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------
#' 2018.08.09.
#' 
#' We are still searching for additional patients in whom to evaluate MRD.
#' 
#' Karen forwarded me the 2017 and 2018 leukemia team patient lists so we
#' can search through them for MRD-positive patients.
#' 
#' First check that no patients are duplicated with the original dataset, 
#' or the joined lists from Olga and Karen.  Then send to Theresa to check
#' consent status and sample availability.
#'-------------------------------------------------------------------------
#'-------------------------------------------------------------------------

require(xlsx); require(dplyr)

#' New patient lists.
new2017 <- read.xlsx('Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/Leukemia team patient list_2017.xlsx',
                     header = TRUE, stringsAsFactors = FALSE, sheetName = '2017', rowIndex = 1:65)
new2018 <- read.xlsx('Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/Leukemia team patient list_2018.xlsx',
                     header = TRUE, stringsAsFactors = FALSE, sheetName = 'all 2018 new patients', rowIndex = 1:23)
new.patients <- rbind(new2017, new2018)

#' Old patient lists.
old.patients <- read.xlsx('Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/olga.and.karen.lists.join.xlsx',
                          header = TRUE, sheetIndex = 1, rowIndex = 1:109, stringsAsFactors = FALSE)
load('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.subjects.expanded.identifiers.rdata')

#' will output any new patients duplicated in previous lists.
#' There are none, so Theresa can review all 86 new patients.  What fun for her.
print(semi_join(new.patients, old.patients, by = c('MRN' = 'mrn')))
print(semi_join(new.patients, metab.clin.expanded.ids, by = 'MRN'))












# Scratch paper -----------------------------------------------------------

olga.dups <- semi_join(olga.list, karen.list, 'mrn')
olga.list <- anti_join(olga.list, karen.list, 'mrn')

#' Create a common set of informative variables.
olga.list$source.list <- 'Olga'
karen.list$source.list <- 'Karen'
karen.list$relapse <- karen.list$cohort

joined.list <- rbind(select(karen.list, mrn, last.name, first.name, source.list, mrd, relapse),
                     select(olga.list, mrn, last.name, first.name, source.list, mrd, relapse))
joined.list$checked.with.rtss <- ifelse(joined.list$source.list == 'Karen', 'Yes', "No")

write.xlsx(joined.list, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Replication set/olga.and.karen.lists.join.xlsx', row.names = FALSE)

rm(list = ls()); gc()
