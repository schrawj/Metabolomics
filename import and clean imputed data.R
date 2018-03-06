#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.05.  
#' Importing the actual metabolomics data into R.
#' Unfortunately, the Metabolon report is not in a format that makes this trivial.
#' I played around with melt and cast but gave up after a while.
#' Far easier to copy the cells of interest in Excel and paste them onto a new sheet by 
#' selecting the transpose option from the paste menu.
#' 
#' 2017.05.11.
#' Received a file from Austin that maps the IDs in Metabolon's spreadsheet to the IDs in
#' my clinical data spreadsheet.  
#' TubeID in "Metabolon_relapse_final.xlsx" = CLIENT IDENTIFIER in Metabolon's data.
#' ID (sheet 1) and Patient# (sheet 2) in Metabolon_relapse_final.xlsx" = ID in the
#' "Metabolomics_subjects" files.  Import "Metabolon_relapse_final.xlsx" and append the TubeID
#' column to the clinical data.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

require(xlsx)
require(dplyr)

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')

subj.metadata <- 
  read.xlsx(file='./Metabolon raw data/Deconstructed Metabolon dataset.xlsx',
            sheetName = 'Sample identifiers long format',colIndex=1:12,header=TRUE,stringsAsFactors = FALSE)

colnames(subj.metadata) <- tolower(colnames(subj.metadata))

lcms <- read.xlsx(file='./Deconstructed Metabolon dataset.xlsx',
                  sheetName='transposed metabolite data', header=TRUE)

chemnames <- c('parent.sample.id', as.character(metab.names$biochemical))
colnames(lcms) <- chemnames

metab <- left_join(subj.metadata, lcms, by='parent.sample.id') 

save(metab,
     file='Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolites.v20170505.1.rdata')


client.id <- read.xlsx(file='./Deconstructed Metabolon dataset.xlsx',
                       sheetName = 'Sample identifiers long format', rowIndex = 1:101, header = TRUE)
names(client.id) <- tolower(names(client.id))

client.id <- arrange(
                      rename(client.id, client.id = client.identifier),
                  client.id)

#' Join the ID variables that will ultimately allow us to pair with the clinical data.
justtheids <- select(client.id, client.id, parent.sample.id)
metab <- left_join(metab, justtheids)
metab <- metab[,c(1,987,2:986)]
save(metab,
     file='Y:/Jeremy Schraw/Metabolomics and Relapse project/Datasets/metabolites.v20170511.1.rdata')

justtheids <- select(id.mappings, client.id, id)
metab <- left_join(metab, justtheids, by='client.id')
metab <- metab[,c(1,2,988,3:987)]
save(metab,
     file='Y:/Jeremy Schraw/Metabolomics and Relapse project/Datasets/metabolites.v20170511.2.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' That takes care of the metabolomics dataset.  As soon as we join the "id" variable
#' from "Metabolon_relapse_list_final" to our clinical dataset, we will be done merging
#' the datasets.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

require(xlsx)
setwd('C:/Users/schraw/Downloads/')

#' Read in the cases.
cases <- read.xlsx('./Metabolon_relapse_list_final.xlsx', sheetIndex=1, header = TRUE)

#' Read in the controls.
controls <- read.xlsx('./Metabolon_relapse_list_final.xlsx', sheetIndex = 2, header = TRUE, rowIndex = 1:61)

#' Row 40 is blank.
controls <- controls[c(1:39,41:60),]

#' DOB, Date of DX and initial WBC count might be useful for checking that joins went smoothly later on, so I'll keep those.
#' But I'm not interested in a lot of these variables.
cases <- cases[,c(1:10,13)]
controls <- controls[,c(1:8,10,17)]
controls$relapse <- 0
controls$relapse <- factor(controls$relapse, levels=0, labels='no')
controls <- controls[,c(1:9,11,10)]

#' in the clinical dataset, "id" is lowercase and a factor.
require(dplyr)
cases <- rename(cases, client.id = TubeID,
                id = ID,
                wbc.count.dx = initial.WBC..x10.3.,
                sex = gender)
cases$id <- as.character(cases$id)


controls <- rename(controls, client.id = TubeID,
                   id = Patient.,
                   diagnosis = Diagnosis,
                   diagnosis.date = Dx.date,
                   wbc.count.dx = Initial.WBC,
                   sex = gender)
controls$id <- as.character(controls$id)
identical(names(cases), names(controls))

id.mappings <- arrange(rbind(cases,controls), by = client.id)

#' Change the id variable in the clinical set to character type.
met$id <- as.character(met$id)

#' Files sould be ready to join.  
#' As a QC measure, let's create a file that has the ids, DOBs and dates of DX from each file, join it, 
#' and see if the dates are concordant.
#' If I left_join met to id, I should also be able to see who is missing from the clinical dataset.
met.dates <- select(met, id, dob, dx_date, initial_wbc)
id.dates <- select(id.mappings, id, DOB, diagnosis.date, wbc.count.dx)

dates <- left_join(id.dates, met.dates, by = 'id')
dates$dx.dates.agree <- ifelse(dates$diagnosis.date == dates$dx_date, TRUE, FALSE)
baddates <- filter(dates, dx.dates.agree == FALSE)
print(baddates)

#' I'm satisfied that the IDs in each file describe the same cases.  However, they are not lining up in all instances.
dates <- left_join(id.dates, met.dates, by = 'id')
no.match <- filter(dates, is.na(dob))
l <- no.match$id
check <- filter(id.mappings, id %in% l)
print(check)

met.unmatched <- left_join(met, id.mappings, by='id')
print(met.unmatched[,30:39])
met.unmatched <- filter(met.unmatched, is.na(client.id))


#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.11.
#' We have 22 rows from the metabolon_relapse_list data file that can't match 
#' with a row in the clinical data.
#' We have 21 rows from the clinical data that can't match with a row in the 
#' Metabolon_relapse_list data.
#' By searching for children using their sex, provider name and DOB in EPIC I was able
#' to ascertain that there are kids in the clinical dataset who are not in the metabolite
#' dataset and vice versa.  Im currently looking up the MRNs for the 
#' kids whose clinical data we do need but don't have in EPIC.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

check<- arrange(check, age.at.dx)
met.unmatched <- arrange(met.unmatched, age_dx)

write.xlsx(met.unmatched,
          file='C:/Users/schraw/Downloads/unmatched kids from clinical dataset.xlsx')
write.xlsx(check,
           file='C:/Users/schraw/Downloads/unmatched kids from metabolon dataset.xlsx')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.12.
#' So, it was relatively straightforward to find MRN numbers for the kids in the 
#' Metabolon dataset who were not in our clinical dataset.  I appended their info
#' to the clinical dataset.
#' Now I want to trim down the clinical dataset to only the individuals who are
#' also in the Metabolon dataset.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

setwd('C:/Users/schraw/Downloads/')
require(xlsx)
require(dplyr)

#' Load in the overpopulated clinical dataset.
met <- read.xlsx(file='./Metabolomics_subjects_5-15-17.xlsx',
                 sheetIndex = 1, header = TRUE, rowIndex = 1:122)

met <- met[,1:30]

save(met, file='Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.subjects.v20170512.1.rdata')

#' Load in the metabolite data.
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolites.v20170511.2.rdata")

id <- select(metab, id)
id$id.in.metab <- TRUE
print(id)

#' Filter out rows in the clinical dataset that do not correspond to rows in the metabolite data.
#' There should be 100 rows in the dataset at the end of this step.
met <- left_join(met, id, by='id')
met <- arrange(
          filter(met, id.in.metab==TRUE),
          id)

#' Yahtzee.  Write it to Excel so we can abstract what we need.
write.xlsx(met,
           file='./Metabolomics_subjects_5-12-17.2.xlsx')

met <- select(met, -id.in.metab)
save(met, file='Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.subjects.20170512.2.rdata')


#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.15.
#' Have mostly completed abstraction on the revised dataset.  Just waiting for Karen to 
#' verify the details of some of the more complex cases.  In the interim, let's load it
#' in and look at whatever summary statistics we can in order to compare it to the 
#' Metabolon deliverables.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

setwd('C:/Users/schraw/Downloads/')
require(xlsx)
require(dplyr)
met.new <- read.xlsx(file='./Metabolomics_subjects_5-15-17.xlsx',
                     sheetIndex = 1, header = TRUE, rowIndex = 1:122, stringsAsFactors = FALSE)
met.new <- met.new[,1:32]

met.new <- left_join(met.new, id)
met.new <- subset(met.new, id.in.metab == TRUE)
with(met.new, table(cohort_new))

#' the data frame 'metab' contains data returned from metabolon.  
#' print some summary statistics for comparison with our new
#' clinical dataset.
with(metab, table(gender))
with(metab, table(group))
with(metab, table(race))

#' summary stats from the updated clinical dataset.
with(met.new, table(gender))
with(met.new, table(cohort_new))

#' race is not coded the same but we can check that the number of NHW
#' and Hispanic individuals match up pretty easily.
nhw <- filter(met.new, race=='W' & ethnicity == 'NH')
hispanic <- filter(met.new, ethnicity == 'H')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.15.
#' Results are encouraging.  Metabolon data lists 59 controls and 41 in the relapse-MRD
#' cohort.  We show 41 in the relapse-MRD and 57 in the 'control' cohort, with two patients
#' awaiting assignment.  Number of NHWs matches exactly.  There are slight discrepancies in 
#' the numbers of Hispanics and gender is 59/41 in Metabolon data and 58/42 in clinical data.
#' Metabolon data has 61 Hispanics; we have 58 with one patient awaiting classification.
#' Investigate discrepancies.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

gender.met.new <- select(met.new, id, last.name, first.name, gender, id.in.metab)
gender.metabolon <- select(metab, id, gender.metab = gender)

#' ID "K" is female (per EPIC record).  Incorrectly classified as male in Metabolon data.
genders <- left_join(gender.met.new, gender.metabolon, by = 'id')
genders <- filter(genders, gender != gender.metab)
print(genders)

print(metab[metab$id=='K',])
metab[11,5] <- 'F'

save(metab,
     file='Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolites v20170515.1.rdata')

with(met.new, table(gender))
with(metab, table(gender))

#' On to ethnicity.
met.new.ethnicity <- select(met.new, id, race, ethnicity)
metabolon.ethnicity <- select(metab, id, race)
ethnicities <- left_join(met.new.ethnicity, metabolon.ethnicity, by = 'id')
print(ethnicities)
mismatched.ethnicities <- ethnicities[c(9,67,78,96),]
print(mismatched.ethnicities)

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.15.
#' There are three instances where ethnicities are mismatched, and one Hispanic Asian
#' which bears looking into.  In all likelihood that individual is actually Hispanic
#' AIAN - American Indian/Alaskan Native.
#' 
#' ID "U" is a high-security record Karen handled.  We'll take her at her word that 
#' this individual is NH.
#' 
#' As far as I can tel ID 1134 may actually be Hispanic Asian.  Asked Karen to verify.
#' 
#' UPDATE 2017.05.25: per Karen, this individual is Hispanic white.
#' This is corrected in version Metabolomics_subjects_5-15-17_kr and above.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolites v20170515.1.rdata")

#' For ID 775, EPIC self reported race-ethnicity form indicates NHW.  
#' Metabolon data is incorrect.  
print(metab[46,1:9])
metab[46,9] <- 'Non-Hispanic White'

#' For ID 1300, EPIC self-report NIH race-ethnicity form confirms NHB.  
#' Metabolon data is again incorrect.
print(metab[metab$id=='1300',])
metab[92,9] <- 'Non-Hispanic Black'
print(metab[92,1:9])

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')

save(metab, file='./metabolites v20170525.1.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.25.
#' Karen has returned the clinical dataset with responses to the flags I gave her.
#' These can be found in Metabolomics_subjects_5-15-17_kr.
#' 
#' I'm going to tidy up that file some by deleting unnecessary columns and then import
#' the data.  
#' 
#' I think we've finally got a production quality dataset to work with. 
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

require(xlsx)
require(dplyr)
require(gmodels)
require(stringr)

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')

metab.clin <- read.xlsx(file='./Clinical raw data/Metabolomics_subjects_5-25-17.xlsx',
                        sheetIndex = 1, rowIndex = 1:122, header = TRUE, stringsAsFactors = FALSE)
names(metab.clin)
metab.clin <- metab.clin[,1:31]
save(metab.clin,
      file='./metabolomics.subjects.v20170525.1.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Now we go through, recoding and computing variables as necessary.
#' But before we get too far, let's cut the dataset down to only those sent to Metabolon.
#' Then, check the demographics against Metabolon's again to identify lingering issues.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolites v20170525.1.rdata")

head(names(metab),25)

metabolon.metadata <- metab[,c(3:7,9)]
head(metabolon.metadata)

#' Compute NIH race-ethnicity from clinical data.
metab.clin$race.ethnicity <- ifelse(metab.clin$ethnicity == 'H', 1,
                                    ifelse(metab.clin$ethnicity == 'NH' & metab.clin$race == 'W', 0,
                                           ifelse(metab.clin$ethnicity == 'NH' & metab.clin$race == 'BLACK', 2, 3)))
metab.clin$race.ethnicity <- factor(metab.clin$race.ethnicity, 
                                    levels = c(0:3), 
                                    labels = c('Non-Hispanic White', 'Hispanic', 'Non-Hispanic Black', 'Other'))
print(metab.clin[,c(29,30,32)])

clinical.metadata <- metab.clin[,c(1,5:7,10,23,28,32)]
head(clinical.metadata)

metadata.check <- left_join(metabolon.metadata, clinical.metadata, by = 'id')
head(metadata.check, 20)

#' Perfect concordance for sex: 42 F, 58 M.
with(metadata.check, CrossTable(gender.x, gender.y))

#' One individual listed as Hispanic in Metabolon set, listed as other in clinical data.
with(metadata.check, CrossTable(race, race.ethnicity))

other <- filter(metab.clin, race.ethnicity == 'Other')
print(other)

hisp <- filter(metab, race == 'Hispanic')
hisp$hisp.in.metab <- TRUE
hisp <- select(hisp, id, hisp.in.metab)

other <- left_join(other, hisp)
print(other)

#' The discordant individual is the high-security record Karen abstracted.  
#' She classified this person as non-Hispanic.
metadata.check$age.dx.floor <- floor(metadata.check$age_dx)
print(metadata.check[,c(2,10,14)])

#' Compare ages.
table(metadata.check$age, metadata.check$age.dx.floor)

#' One subject, ID 'E' does not align.
#' Clinical data shows 2005 DOB and tx ongoing; Metabolon data incorrect.
metadata.check$age.diff <- with(metadata.check, age - age.dx.floor)
age.diff <- filter(metadata.check, age.diff != 0)
print(age.diff)

print(metab[metab$id=='E',])
metab[5,4] <- 10

save(metab, file='./metabolites.v20170525.2.rdata')

#' Compare relapse-MRD classifications.
#' 2 kids are discordant.
with(metadata.check, CrossTable(group, cohort_new))

#' ID 'J' is listed as relapsed in the Metabolon set but not the clinical set.
#' Karen reviewed this one for cytogenetics and didn't correct me on relapse.
#' TODO: triple check this when you get EPIC access again.
classifications <- filter(metadata.check, cohort_new == 'no relapse' & group == 'Case')
print(classifications)

#' ID '1125' is listed as relapsed in the clinical dataset but not the Metabolon dataset.
#' This child relapsed after the initial data collection period.  Relapse was discovered 
#' when I went back through EPIC in 4/2017.
classifications <- filter(metadata.check, cohort_new == 'relapse MRD' & group == 'Control')
print(classifications)

print(metab[metab$id=='1125',])
metab[80,6] <- 'Case'

save(metab, file='./metabolites.v20170525.2.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' That's about all the metadata we can compare.  Given occasional errors in the 
#' Metabolon metadata, we'll use the clinical data instead. 
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

str(metab.clin)

#' There are a lot of identifiers in this dataset.  More than we need or should have.
#' Save them in  another dataset in case we need them then delete them.
metab.clin.expanded.ids <- metab.clin[,c(1:5)]
save(metab.clin.expanded.ids, file='./metabolomics.subjects.expanded.identifier.list.v20170525.1.rdata')

#' Remove extraneous identifiers.
metab.clin <- select(metab.clin, -RTSS.ID, -last.name, -first.name, -MRN)

#' Convert 'diagnosis' to a factor named 'immunophenotype'.
metab.clin$diagnosis <- substr(metab.clin$diagnosis, 1, 5)
metab.clin$immunophenotype <- with(metab.clin, ifelse(diagnosis == 'B-ALL',0,
                                                      ifelse(diagnosis == 'T-ALL',1,2)))
metab.clin$immunophenotype <- factor(metab.clin$immunophenotype, 
                                     levels = c(0:2),
                                     labels = c('B-lineage','T-lineage','This must be an error'))
table(metab.clin$diagnosis, metab.clin$immunophenotype)

metab.clin <- select(metab.clin, -diagnosis)
metab.clin <- metab.clin[,c(1,2,28,3:27)]

#' The cytogenetics variable still lists all the multiple rearrangements for each case.
#' Trim this down to the primary aberration and compute the 'neutral', 'favorable' and 
#' 'unfavorable' groups we discussed previously.
metab.clin.expanded.cytogenetics <- select(metab.clin, id, cytogenetics)

save(metab.clin.expanded.cytogenetics, 
     file='./metabolomics.subjects.expanded.cytogenetics.v20170525.1.rdata')

#' Note that this code keeps text to the left of semicolons.
#' There are semicolons in the karyotype data, but it doesn't matted because
#' the prognostic/recurrent alterations are referred to by their gene symbols
#' instead of their chromosomal coordinates.
metab.clin$cyto.simple <- gsub(';.*$','',metab.clin$cytogenetics)

#' Per the 2017.05.12 meeting notes we decided we should only have 'favorable/neutral'
#' and 'unfavorable' groups.  For the moment I will compute a 3-group classifier and 
#' then compress it to 2 groups.
metab.clin$cyto.group.threecat <- with(metab.clin, ifelse(cyto.simple == 'ETV6-RUNX1' | cyto.simple == 'hyperdiploid', 1,
                                                          ifelse(cyto.simple == 'BCR-ABL' | cyto.simple == 'MLL rearrangement' | cyto.simple == 'iAMP21', 2,
                                                                 ifelse(cyto.simple == 'not performed', 3, 0))))
metab.clin$cyto.group.threecat <- factor(metab.clin$cyto.group.threecat, 
                                         levels=c(0:3),
                                         labels=c('Neutral cytogenetics', 'Favorable cytogenetics', 
                                                  'Unfavorable cytogenetics', 'Cytogenetics not performed'))
print(metab.clin[,c(1,29,30)])

metab.clin$cyto.group.twocat <- with(metab.clin, ifelse(cyto.group.threecat == 'Neutral cytogenetics' | cyto.group.threecat == 'Favorable cytogenetics', 0,
                                                        ifelse(cyto.group.threecat == 'Unfavorable cytogenetics', 1, 2)))
metab.clin$cyto.group.twocat <- factor(metab.clin$cyto.group.twocat, 
                                       levels=c(0:2),
                                       labels=c('Favorable or neutral', 'Unfavorable', 'Cytogenetics not performed'))
print(metab.clin[,c(1,29:31)])

metab.clin <- select(metab.clin, -cytogenetics)
metab.clin <- metab.clin[,c(1,23,2,24,25,27,3,28:30,4:22,26)]

#' NCI risk group needs to be computed for the new additions to the dataset.
#' The NCI risk catgories were written for kids with B-cell ALL, but we will compute it for all children, including the 12
#' with T-lineage ALL.
metab.clin$nci.risk.group <- with(metab.clin, ifelse(initial_wbc >= 50 | age_dx >= 10, 1,
                                                     ifelse(initial_wbc < 50 & age_dx < 10, 0, 2)))
metab.clin$nci.risk.group <- factor(metab.clin$nci.risk.group,
                                    levels = c(0:2),
                                    labels = c('Standard risk', 'High risk', 'This is an error'))

metab.clin <- select(metab.clin, -nci_risk_group, -nci_risk_new)
metab.clin <- metab.clin[,c(1:13,29,14:28)]

identical(metab.clin$cohort, metab.clin$cohort_new)

#' 'protocol' is probably another variable that can safely be moved to our identifiers dataset.
protocol <- select(metab.clin, id, protocol)
metab.clin.expanded.ids <- left_join(metab.clin.expanded.ids, protocol, by = 'id')
save(metab.clin.expanded.ids,
     file='metabolomics.subjects.expanded.identifiers.v20170525.1.rdata')

metab.clin <- select(metab.clin, -protocol)

save(metab.clin, 
     file='./metabolomics.subjects.v20170525.2.rdata')

#' 'fu_time' appears to be a string of NAs.  Recalculate.
metab.clin$fu.time <- as.numeric(metab.clin$last_fu_date - metab.clin$dx_date)
summary(metab.clin$fu.time)

metab.clin <- select(metab.clin, - fu_time)
metab.clin <- metab.clin[,c(1:11,20,28,12:19,21:27)]

#' Convert gender to a factor named sex.
metab.clin$sex <- with(metab.clin, ifelse(gender=='M', 2, 1))
metab.clin$sex <- factor(metab.clin$sex, levels=c(1, 2), labels=c('Female','Male'))

metab.clin <- select(metab.clin, -gender)
metab.clin <- metab.clin[,c(1,28,2:27)]

save(metab.clin, 
     file='./metabolomics.subjects.v20170525.3.rdata')


#' dob is problematic.  The majority of cells are the weird Excel 'number of days since antiquity' format.
#' A few are dates, but with two digit years instead of four.
#' This is a pain in the dick.  I think I have to pull out all the rows with weird dates, operate on them 
#' individually and bind them back in.  Like some kind of peasant.
metab.clin <- arrange(metab.clin, dob)
print(metab.clin[,c(1,3)])

datefix <- metab.clin[c(16:89,94:104),]
datefix$dob <- as.Date(as.numeric(datefix$dob), origin='1899-12-30')
print(datefix[,c(1,3)])

#' Choose five and compare them to the Excel sheet.  Rows 6,11,29,51,83.
#' All five match exactly.

#' Pick out the kids with two digit years.
datefix2 <- metab.clin[c(20,21,23,25,33),]
birthdates <- c('2011-05-24','2006-05-27','2011-06-22','2011-06-08','1999-09-19')
datefix2$dob <- birthdates
datefix2$dob <- as.Date(datefix2$dob)


#' and finally the kids whose dob is a character formatted as a date.
datefix3 <- metab.clin[c(1:19,22,24,26:32,34:36),]
datefix3$dob <- as.Date(datefix3$dob, format = '%m/%d/%Y')
datefix3[14,3] <- as.Date('1998-03-12')

metab.clin <- rbind(datefix, rbind(datefix2, datefix3))

#' On my way out of the office for evening, but rename variables first.  No underscores.
metab.clin <- rename(metab.clin,
                     dx.date = dx_date,
                     last.fu.date = last_fu_date,
                     age.dx = age_dx,
                     init.wbc = initial_wbc,
                     eotdate.notes = eotdate_notes,
                     day29mrd.notes = day29mrd_notes,
                     was.flagged = flag_for_Karen,
                     relapse.type = relapse_type,
                     date.of.relapse = Date_of_relapse,
                     relapse.type.notes = relapse_type_notes)

save(metab.clin, file='./metabolomics.subjects.v20170525.4.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.26.
#' Change 'relapse' to a factor: yes, no, < 3 years follow up.
#' Break subjects into two groups: 1) alive/dead and 2) alive, treatment ongoing.
#' For kids on therapy with < 3 years follow-up, set relapse to NA.
#' There is one kid with ~1060 days of follow up.  
#' Included them.  Manually computed their relapse value.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

done.w.tx <- filter(metab.clin, alive.dead == 'alive' | alive.dead == 'dead')
still.in.tx <- filter(metab.clin, alive.dead == 'alive; tx ongoing')

#' Computing relapse is more straightforward in the kids who are off therapy.  First clean up NAs.
which(is.na(done.w.tx$relapse))
#' These are the two kids with second primary tumors.  Karen confirmed neither is in our relapse cohort.
print(done.w.tx[c(89,94),])
done.w.tx[c(89,94),17] <- 'no'

#' Convert to a factor.
done.w.tx$relapse <- with(done.w.tx, ifelse(relapse=='no', 0, 1))

done.w.tx$relapse <- factor(done.w.tx$relapse, 
                            levels = c(0, 1), labels = c('No Relapse', "Relapse"))

#' For kids still on therapy, only compute an outcome for relapse is there is at least 3 years follow-up.
still.in.tx$relapse <- with(still.in.tx, ifelse(fu.time >= 1095 & relapse == 'no', 0,
                                                ifelse(fu.time >= 1095 & relapse == 'yes', 1, NA)))
still.in.tx$relapse <- factor(still.in.tx$relapse, 
                              levels = c(0, 1),
                              labels = c('No Relapse', "Relapse"))

#' The child with borderline follow-up time.
still.in.tx[9,17] <- 'No Relapse'

metab.clin <- rbind(done.w.tx, still.in.tx)

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' There's an error in my relapse variable.  ID 'E' relapsed at < 3 years and was 
#' incorrectly filed as NA.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

metab.clin[metab.clin$id=='E',17] <- 'Relapse'

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' The MRD variable is messy.  Some rows have text values.
#' The rows with numerica values were divided by 100 as R converted them from a
#' percentage to numeric.
#' 
#' I'll apply the same logic as for DOB: split into groups, one with text values and one
#' with numeric values, then join back together.
#' 
#' Our endpoint is to have MRD defined as a factor with levels 'yes' corresponding to 
#' day 29 MRD by flow cytometry >= 0.01% or positive findings for FISH or positive 
#' findings for karyotype and 'no' corresponding to day 29 MRD < 0.01% by flow.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

#' The lay of the land. 
metab.clin <- arrange(metab.clin, id)
print(metab.clin[,c(1,22)])

#' Subset the data.
mrd.numeric <- subset(metab.clin[c(1:59,63:70,73,75:82,84:88,90,91,93:121),])
mrd.numeric$mrd <- as.numeric(mrd.numeric$day29mrd)
mrd.numeric$mrd <- mrd.numeric$mrd*100
mrd.numeric$mrd.cat <- factor(ifelse(mrd.numeric$mrd < 0.01, 0, 1),
                            levels = c(0, 1),
                            labels = c('MRD negative', 'MRD positive'))
mrd.numeric <- select(mrd.numeric, -mrd)
                            
l <- c(mrd.numeric$id)
mrd.text <- subset(metab.clin, !(metab.clin$id %in% l))
print(mrd.text[,c(1,22)])

#' As it happens the only two scenarios in this group are that MRD was not assessed, or 
#' the patient was below the cutoff but not quantified exactly.
mrd.text$mrd.cat <- factor(
  
                            ifelse(mrd.text$day29mrd == 'MRD not done', NA, 0),
                            
                        levels = c(0, 1),
                        labels = c('MRD negative', 'MRD positive'))

metab.clin <- rbind(mrd.numeric, mrd.text)

save(metab.clin,
     file = './metabolomics.subjects.v20170526.1.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' I inherited a dataset with multiple columns named 'cohort'.
#' I may need to compute a completely new one anyways now that I've computed my own
#' versions of the MRD and relapse variables.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

metab.clin$group.fine <- with(metab.clin,
                          
                          factor(
  
                              ifelse(relapse == 'No Relapse' & mrd.cat == 'MRD negative', 0,
                                     ifelse(relapse == 'No Relapse' & mrd.cat == 'MRD positive', 1,
                                            ifelse(relapse == 'Relapse' & mrd.cat == 'MRD negative', 2,
                                                   ifelse(relapse == 'Relapse' & mrd.cat == 'MRD positive', 3, NA)))),
                              
                          levels = c(0:3),
                          labels = c('MRD- Relapse-', 'MRD+ Relapse-', 'MRD- Relapse+', 'MRD+ Relapse+')))

metab.clin$group.coarse <- with(metab.clin,
                                
                                factor(
                                  
                                    ifelse(relapse == 'Relapse' | mrd.cat == 'MRD positive', 1,
                                           ifelse(relapse == 'No Relapse' & mrd.cat == 'MRD negative', 0, NA)),
                                
                                levels = c(0, 1),
                                labels = c('No relapse or MRD', 'Relapse-MRD')))

metab.clin <- select(metab.clin, -cohort, -cohort_new)
metab.clin <- metab.clin[,c(1,3,2,4:6,15,16,14,7:13,21,22,27,17,23:25,28,29,18:20,26)]
metab.clin <- metab.clin[,c(1:6,9,7,8,10:29)]

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Convert 'alive.dead' to a factor.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

unique(metab.clin$alive.dead)

metab.clin$alive.dead <- with(metab.clin, 
                              factor(
                                  ifelse(alive.dead == 'dead', 2,
                                         ifelse(alive.dead == 'alive; tx ongoing', 1, 0)),
                              
                              levels = c(0:2),
                              labels = c('Alive', "Alive; on therapy", 'Deceased')))

save(metab.clin, 
     file='./metabolomics.subjects.v20170526.2.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Relapse type needs to be cleaned.  
#' 1) Convert to a factor.
#' 2) Account for length of follow up.  
#' 3) There is one child who I show did not relapse but for whom relapse type is NA.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

print(metab.clin[,c(1,20,22,23)])

#' Fix problem 3.
metab.clin[metab.clin$id=='U',22] <- 'no relapse'

#' The six kids with 'relapse' == NA and relapse.type == no.relapse are the ones with
#' inadequate follow up.
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.subjects.expanded.identifiers.v20170525.1.rdata")

l <-c('B','C','D','F','G','H')

followup <- filter(metab.clin, id %in% l)
followup$relapse.type <- NA

metab.clin <- filter(metab.clin, !(id %in% l))

metab.clin <- rbind(metab.clin, followup)

#' Convert to a factor.
metab.clin$relapse.type <- with(metab.clin,
                                
                                factor(
                                  
                                  ifelse(relapse.type == 'no relapse', 0,
                                         ifelse(relapse.type == 'isolated BM', 1,
                                                ifelse(relapse.type == 'isolated extramedullary', 2,
                                                       ifelse(relapse.type == 'combined', 3, NA)))),
                                  
                                levels = c(0:3),
                                labels = c('No relapse', 'Isolated BM',
                                           'Isolated extramedullary', 'Combined')))

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' There are still some clinical variables that would be good to have for reference,
#' but which I don't anticipate using in the analysis. 
#' Move them to a separate data frame.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

metab.clin.expanded.clinical.outcomes <- select(metab.clin, id, day29mrd, day29mrd.notes, date.of.relapse, 
                                                            relapse.type, relapse.type.notes, eotdate, eotdate.notes,
                                                            last.fu.date, was.flagged)

save(metab.clin.expanded.clinical.outcomes,
     file='./metabolomics.subjects.expanded.clinical.outcomes.rdata')

#' Compute a time to relapse variable and then get all these dates out of the dataset.
print(metab.clin[,c(1,14,21)])

metab.clin$time.to.relapse <- as.numeric(metab.clin$date.of.relapse - metab.clin$dx.date)

print(metab.clin[,c(1,14,21,31)])

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Some final data cleaning and sorting.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

metab.clin <- select(metab.clin, -day29mrd.notes, -date.of.relapse, -relapse.type.notes, -eotdate, -eotdate.notes, 
                     -was.flagged, -last.fu.date)

metab.clin <- metab.clin[,c(1:19,23,20:22)]

metab.clin <- rename(metab.clin,
                     race.eth = race.ethnicity,
                     nci.risk = nci.risk.group,
                     cyto = cyto.simple,
                     cyto.three.cat = cyto.group.threecat,
                     cyto.two.cat = cyto.group.twocat,
                     mrd = mrd.cat)

metab.clin <- select(metab.clin, -day29mrd)

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' I might as well repair the truncated cytogenetics values for kids with less common
#' translocations.  It's the one black eye on my otherwise impeccable dataset.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

metab.clin[metab.clin$id=='1136',11] <- 't(1;3)'
metab.clin[metab.clin$id=='Y',11] <- 't(1;2)'
metab.clin[metab.clin$id=='H',11] <- 't(13;16)'
metab.clin[metab.clin$id=='907',11] <- 'possible t(9;9)'
print(metab.clin[,c(1,11)])

metab.clin <- arrange(metab.clin, id)

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Lastly, we should pare it down to the 100 kids with metabolite data.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

load("//discovery1/bcm-dldcc-epi$/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolites.v20170525.2.rdata")

id <- select(metab, id)

metab.clin <- left_join(id, metab.clin, by = 'id')

save(metab.clin, 
     file='./metabolomics.subjects.v20170526.3.rdata')

write.xlsx(metab.clin,
           file='./metabolomics.subjects.v20170526.xlsx')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' Metabolon metadata is occasionally wrong.  Remove it, and join to the clinical file.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

metab <- select(metab, -parent.sample.id, -client.id, -age, -gender, -group, -group.number, -race)

metab.bio <- metab

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')

save(metab.bio,
     file='./metabolites.v20170526.1.rdata')

met <- left_join(metab.clin, metab.bio, by = 'id')

save(met,
     file='metabolomics.relapse.v20170526.1.rdata')

write.xlsx(met,
           file='./metabolomics.relapse.v20170526.xlsx')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.30.
#' I performed a large number of operations on the clinical data.  It wouldn't hurt to 
#' compare the metadata one last time to make sure noting horrible has happened.
#' That said, it can wait until after the weekend.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/metabolomics.subjects.v20170526.3.rdata')
#' Have to load in an older Metabolon dataset.  Final one has the metadata scrubbed.
load('./Datasets/Old R datasets/metabolites.v20170525.2.rdata')

require(dplyr)

clinical.metadata <- select(metab.clin, id, age.dx, sex, relapse, mrd, group.coarse, group.fine, race, ethnicity, race.eth)
metab.metadata <- select(metab, parent.sample.id, client.id, id, age, gender, group, group.number, race)
metadata.comparison <- left_join(clinical.metadata, metab.metadata, by = 'id')

require(gmodels)

with(metadata.comparison, CrossTable(sex, gender))

summary(metadata.comparison$age)
summary(metadata.comparison$age.dx)

#' One child Hispanic in Metabolon set, 'Other' in clincal set.
with(metadata.comparison, CrossTable(race.eth, race.y))

race.eth <- left_join(
  
                        select(metab.clin, id, race.eth),
                        select(metab, id, race),
                      
                      by = 'id')

race.eth$same <- as.character(race.eth$race.eth) == as.character(race.eth$race)

#' That child is ID 'U', whom Karen handled.  She says NH.
print(race.eth)

#' Fix the Metabolon data.
metab[metab$id=='U',9] <- 'Non-Hispanic Asian'

#' This is a problem.
with(metadata.comparison, CrossTable(group.coarse, group), na.rm = TRUE)

groups <- arrange(
  
                  left_join(
        
                    select(metab.clin, id, relapse, mrd, group.coarse),
                    select(metadata.comparison, id, group, group.number),
                    
                  by = 'id'),
              
          group)
            
print(groups)

#' ID 'J' is no relapse & MRD negative in the clinical set; a case in the Metabolon data.
#' clinical spreadsheet confrims no relapse & MRD negative.

metab[metab$id=='J', 6] <- 'Control'
metab[metab$id=='J', 1:6]

#' There are six kids classified as relapse-MRD in the clinical data and controls in the Metabolon data.
#' IDs: 
#'    923  - 0.02% MRD. Same as in Metabolon_relapse_final.   
#'    1004 - 3% MRD. Given as 0.0% in Metabolon_relapse_final.  
#'    1025 - 0.3% MRD. Same as in Metabolon_relapse_final.
#'    1089 - 0.8% MRD. Same as in Metabolon_relapse_final.  
#'    1169 - 0.05% MRD. Same as in Metabolon_relapse_final.
#'    1251 - 1.3% MRD. Same as in Metabolon_relapse_final.
#'    
#' So, it doesn't look like the issue is with the integrity of the data, except maybe for ID 1004.
#' Without knowing the speciifics my assumption is that there was originally a different definition of 
#' 'case' and 'control' than the one we agreed upon. 

problem.child <- subset(metadata.comparison, group.coarse == 'Relapse-MRD' & group == 'Control')
print(problem.child)

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' A few final cosmetic changes.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

metab.clin <- rename(metab.clin,
                     eth = ethnicity)

save(metab.clin, 
     file='./Datasets/metabolomics.subjects.v20170530.1.rdata')

load("./Datasets/metabolomics.relapse.v20170526.1.rdata")

met <- rename(met, 
              eth = ethnicity)

save(met,
     file = './Datasets/metabolomics.relapse.v20170530.1.rdata')

#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.05.31.
#'
#' One thing I did not check is whether the cells housing the metabolite data have values 
#' are correct.  I.e., does the cell for subject X & metabolite Y actually correspond to
#' the value for that same combination in the Metabolon raw data?
#' 
#' I noticed the metbolites are listed in a different order in this dataset than in the 
#' raw data and am concerned as to why that might be.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------

require(dplyr)
require(xlsx)
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.subjects.expanded.identifiers.rdata")
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20170531.1.rdata")
setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')

check <- left_join(met, 
                   
                   select(metab.clin.expanded.ids, id, client.identifier, by = 'id'))

check <- check[,c(1,1004,2:1003)]

#' Metabolite measurements begin at column 27.
head(names(check),30)

check[1:20,c(2,27:38)]

#' FUUUUUUUUUCKK.  They don't match the raw data.  
#' It looks like something went wrong when I named the columns in met.
#' Let's try again.
names <- read.xlsx(file='./Datasets/Metabolon raw data/Original Metabolon dataset.xlsx',
                   sheetName = 'ScaledImpData', colIndex = 2, rowIndex = 13:990, header = TRUE,
                   stringsAsFactors = FALSE)

names <- c(names$BIOCHEMICAL)

metabolites <- met[,27:1003] 
colnames(metabolites) <- names

clin <- met[,1:26]

met <- cbind(clin, metabolites)

#' I think that should have repaired it.  Let's check the values again.
check <- left_join(met, 
                   
                   select(metab.clin.expanded.ids, id, client.identifier),
                   
              by = 'id')

check <- check[,c(1,1004,2:1003)]

#' Yep, that's better.
save(met,
     file='./Datasets/metabolomics.relapse.v20170531.2.rdata')

write.xlsx(met,
           file='./Datasets/metabolomics.relapse.v20170531.2.xlsx')

require(ggplot2)

#' Box plots confirm this issue is fixed.
gg <- ggplot(data = met) + geom_boxplot(aes(x=group.coarse, y=met$"arachidonate (20:4n6)"))
print(gg)

gg <- ggplot(data = met) + geom_boxplot(aes(x=group.coarse, y=met$"arginine"))
print(gg)



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.06.19.
#' 
#' We sent a list of the high scoring unnamed metabolites to Metabolon and had one come 
#' back named: x-23665 has been identified as N-carboxyethylphenylalanine in the time 
#' since these samples were analyzed.  Update its name and annotation.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
require(dplyr)
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20170616.1.rdata")
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.candidate.compounds.v20170616.1.rdata")
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.compounds.expanded.info.v20170616.rdata")
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/p.values.and.mean.ratios.split.method.v20170616.rdata")
setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')

met <- rename(met, `N-carboxyethylphenylalanine` = `X - 23665`)
candidates <- rename(candidates, `N-carboxyethylphenylalanine` = `X - 23665`)
metabolite.annotations[914,1] <- 'N-carboxyethylphenylalanine'
metabolite.annotations[914,2] <- 'Amino Acid'
p.values.split.applied[748,1] <- 'N-carboxyethylphenylalanine'
p.values.split.applied[748,8] <- 'Amino Acid'

save(met, 
     file = './metabolomics.relapse.v20170619.1.rdata')

save(candidates, 
     file = './metabolomics.relapse.candidate.compounds.v20170619.1.rdata')

save(metabolite.annotations, 
     file = './Expanded datasets/metabolomics.compounds.expanded.info.v20170619.rdata')

save(p.values.split.applied,
     file = './Expanded datasets/p.values.and.mean.ratios.split.method.v20170619.rdata')



#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
#' 2017.07.18.
#' 
#' At meeting yesterday with Michael, Philip and Karen we decided to filter the unnamed
#' compounds.
#'---------------------------------------------------------------------------------------
#'---------------------------------------------------------------------------------------
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20170619.1.rdata")
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.compounds.original.scale.v20170619.1.rdata")
load("Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/Expanded datasets/metabolomics.compounds.expanded.info.v20170716.rdata")
setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')

tmp <- filter(metabolite.annotations, is.na(super.pathway))
tmp <- c(tmp$compound)
met <- met[, !(colnames(met) %in% tmp)]
met.orig <- met.orig[, !(colnames(met.orig) %in% tmp)]

#' Neither data frame contains any columns with text patterns matching the unknown compounds.
tmp <- grep('X - ', colnames(met))
tmp <- grep('X - ', colnames(met.orig))

save(met, file = './metabolomics.relapse.v20170718.rdata')
save(met.orig, file = './metabolomics.relapse.original.scale.v20170718.rdata')




# Update an incorrect DOB -------------------------------------------------

met[met$id == 'O','dob'] <- as.Date('2001-06-08')
met[met$id == 'O','age.dx'] <- 12.35342

save(met, file = 'Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/metabolomics.relapse.v20180214.1.rdata')




# Update an incorrect weight ----------------------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')
load('metabolomics.relapse.v20180216.1.rdata')
load('./Expanded datasets/metabolomics.subjects.expanded.anthropometrics.rdata')

#' Weight for one subject was identified as being high by 1 kg.
met[met$id == '1396', 'weight.kg.dx'] <- 25.2

anthro <- select(anthro, -height.cm.dx.y)
anthro <- rename(anthro, age.dx = age.dx.x, height.cm.dx = height.cm.dx.x, age.dx = age.dx.y)
anthro[anthro$id == '1396', 'bmiz'] <- 1.77
anthro[anthro$id == '1396', 'bmipct'] <- 96

save(met, file = 'metabolomics.relapse.v20180219.1.rdata')
save(anthro, file = './Expanded datasets/metabolomics.subjects.expanded.anthropometrics.rdata')




# Compute a categorical variable for BMI ----------------------------------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/Datasets/')
load('metabolomics.relapse.v20180219.1.rdata')

#' Using CDC cutoffs:
#' https://www.cdc.gov/healthyweight/assessing/bmi/childrens_bmi/about_childrens_bmi.html
met$weight.status <- factor(
                            ifelse(met$bmipct < 5, 1,
                            ifelse(met$bmipct >= 5 & met$bmipct < 85, 0,
                            ifelse(met$bmipct >= 85 & met$bmipct < 95, 2,
                            ifelse(met$bmipct >= 95, 3, NA)))),
                            
                            levels = c(0:3),
                            labels = c('Normal Weight','Underweight','Overweight','Obese'))

met$overwt.status <- factor(
                            ifelse(met$bmipct < 85, 0,
                            ifelse(met$bmipct >= 85, 1, NA)),
                            
                            levels = c(0,1),
                            labels = c('Not Overweight/Obese', 'Overweight/Obese'))

met <- met[,c(1:13,745,746,14:744)]

save(met, file = 'metabolomics.relapse.v20180222.1.rdata')



# Three additional children met follow-up requirement for relapse ---------

setwd('Y:/Jeremy Schraw/Metabolomics and relapse project/')
load('./Datasets/metabolomics.relapse.v20180222.1.rdata')
load('./Datasets/Expanded datasets/metabolomics.subjects.expanded.identifiers.rdata')

tmp <- met[is.na(met$relapse), ]
print(tmp[, c('id','dx.date')]) # IDs are B, C, D.

tmp <- filter(metab.clin.expanded.ids, id %in% c('B','C','D'))

#' Upon review of most recent EPIC data (as of 03/06/2018), all three are still in remission.
#' All three are still on therapy.

ids <- c('B','C','D')

tmp <- filter(met, id %in% ids)
met <- filter(met, !(id %in% ids))

for (i in tmp){
  tmp$relapse <- as.factor('No Relapse')
  tmp$fu.time <- as.Date('2018-03-06') - tmp$dx.date
}

met <- rbind(met, tmp)

save(met, file = './Datasets/metabolomics.relapse.v20180306.1.rdata')
