# Name: 03_downloadData_GSE83147.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 23/10/2017
# Desc: Download the relevant datasets, normalise and create database entries

source('header.R')

library(GEOquery)
library(downloader)

## open the soft format and raw data
gse =  getGEO(filename = 'dataExternal/gse83147/GSE83147_series_matrix.txt.gz')

# get the samples from the expression set object
dfSamples = pData(gse)

# create database entries
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and metafile table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe MetaFile;'))
cFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

## create the entry for samples
cSampleCol
g_pid = 9
g_did = 18
title = paste(as.character(dfSamples$title), rownames(dfSamples))
description = paste(as.character(dfSamples$source_name_ch1), as.character(dfSamples$growth_protocol_ch1), sep='. ')
f = grepl('trauma', as.character(dfSamples$characteristics_ch1))
group1 = rep('RA', times=6)
group1[f] = 'Trauma'
group2 = gsub('Sex: (\\w+)', '\\1', as.character(dfSamples$characteristics_ch1.1))
group3 = gsub('\\D', '', as.character(dfSamples$characteristics_ch1.2))

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=title, description=description, group1=group1, 
                       group2=group2, group3=group3)
# write this data to the database
rownames(dfSamples) = NULL

# # write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)

## check annotation data
str(fData(gse))
## NOTE: no annotation available
## store the grouping factors in affy object
gse$fCondition = dfSamples$group1
gse$fSerum = dfSamples$group2
gse$fAge = as.numeric(as.character(dfSamples$group3))

## check normalisation 
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

m = exprs(gse)
range(m)
#m = m + abs(min(m)) + 0.1
oDiag.1 = CDiagnosticPlots(m, 'Soft format')

fBatch = gse$fCondition

## check normalisation
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

## this should give an error as scaling can't be done
## if all the vector 0 for PCA
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)

plot.PCA(oDiag.1.2, fBatch, legend.pos = 'topright')

plot.dendogram(oDiag.1.2, fBatch, labels_cex = 0.8, cex.main=0.7)

### save the object and create metafile entry
getwd()
n = make.names(paste('GSE83147 fibroblast-like synoviocytes from rheumatoid arthritis patients and trauma patients.rds'))
n2 = paste0('~/Data/MetaData/', n)
save(gse, file=n2)

## access the mysql database to create entry for this data object
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='GSE83147 fibroblast-like synoviocytes from rheumatoid arthritis patients and trauma patients. Data is normalised but no annotation found so far and includes non-coding rna as well.')
#dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

