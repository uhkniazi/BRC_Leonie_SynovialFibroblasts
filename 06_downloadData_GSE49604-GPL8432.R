# Name: 06_downloadData_GSE49604-GPL8432.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 24/10/2017
# Desc: Download the relevant datasets, normalise and create database entries

source('header.R')

library(GEOquery)
library(downloader)

## open the soft format and raw data
oExp =  getGEO(filename = 'dataExternal/gse49604/GSE49604-GPL8432_series_matrix.txt.gz')

## read raw data files and normalize
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(lumi)

# add lumi nuIDs 
oExp = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )
dim(oExp)

i = rowSums(exprs(oExp))
table(is.na(i))

## data normalization
# normalize and log2 transform the data using lumi
# remove negative values first and set minimum value to 1
range(exprs(oExp))

# oExp = lumiN(oExp, method='rsn')
# gc()

# get the samples from the expression set object
dfSamples = pData(oExp)

# sanity check
identical(colnames(oExp), rownames(pData(oExp)))

# create database entries
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and metafile table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

## create the entry for samples
cSampleCol
g_pid = 9
g_did = 23
title = paste(as.character(dfSamples$title), rownames(dfSamples))
description = paste(as.character(dfSamples$source_name_ch1), as.character(dfSamples$characteristics_ch1),
                    as.character(dfSamples$characteristics_ch1.2),
                    as.character(dfSamples$characteristics_ch1.3),
                    as.character(dfSamples$characteristics_ch1))
group1 = rep('RA', times=nrow(dfSamples))
f = grepl('Healthy', dfSamples$source_name_ch1)
group1[f] = 'Healthy'
group2 = rep('Blood-Monocytes', times=nrow(dfSamples))
f = grepl('RA', dfSamples$source_name_ch1)
group2[f] = 'Synovial-Macrophages'

group3 = rep('female', times=nrow(dfSamples))

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=title, description=description, group1=group1, 
                       group2=group2, group3=group3)
# write this data to the database
rownames(dfSamples) = NULL

# # write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)

## annotation data
str(fData(oExp))
## store the grouping factors in affy object
oExp$fCondition = dfSamples$group1
oExp$fGender = dfSamples$group3
oExp$fTissue = dfSamples$group2

## check normalisation 
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(exprs(oExp), 'Soft format')

fBatch = oExp$fCondition

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

### save the x.affy object and create metafile entry
getwd()
n = make.names(paste('GSE49604-GPL8432 Healthy PBMC and RA macrophages.rds'))
n2 = paste0('~/Data/MetaData/', n)
save(oExp, file=n2)

## access the mysql database to create entry for this data object
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='GSE49604-GPL8432 Healthy PBMC and RA macrophages. This study has another data set check for data id 21.')
#dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

