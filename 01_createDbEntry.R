# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 25/07/2018


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe MetaFile;'))
cFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

## load the metadata file
dfMeta = read.csv('dataExternal/ProjectID_15_sampleSheet.csv', header=T, stringsAsFactors = F)
str(dfMeta)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta$Title, 
                       description= dfMeta$Description,
                       group1 = dfMeta$Grouping, group2= dfMeta$Adjust1)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 30;'))

cFileCol

dfCounts = read.csv('dataExternal/ProjectID_15_rawData.csv', header = T, sep=',', stringsAsFactors = F)
dim(dfCounts)
str(dfCounts)
cn = colnames(dfCounts)
cn = cn[-c(1:8)]
t = dfSamples$title

## match the names and ordering in sample info and raw data
cn = gsub('X(\\d+).+', replacement = '\\1', cn)
cn = gsub('.Recovered.', '', cn)

t = gsub('-.+', '', t)
t = gsub('\\(Recovered\\)', '', t)

identical(cn, t)
# # put the names in the same order
# table(colnames(dfCounts) %in% dfSamples$title)
# 
# i = match(dfSamples$title, colnames(dfCounts))
# 
# identical(dfSamples$title, colnames(dfCounts)[i])
# dfCounts = dfCounts[,i]

# # sanity check
# identical(dfSamples$title, colnames(dfCounts))

n = make.names(paste('dfCounts for Proteomics data Jamie data id 30 rds'))
n2 = paste0('~/Data/MetaData/', n)
save(dfCounts, file=n2)

dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='data frame of proteomics count data for Jamie data id 30 rds')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)
