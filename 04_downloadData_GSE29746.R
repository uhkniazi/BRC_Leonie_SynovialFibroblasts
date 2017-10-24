# Name: 04_downloadData_GSE29746.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 23/10/2017
# Desc: Download the relevant datasets, normalise and create database entries

source('header.R')

library(GEOquery)
library(downloader)
library(limma)
## open the soft format and raw data
gse =  getGEO(filename = 'dataExternal/gse29746/GSE29746_series_matrix.txt.gz')

## read raw data files and normalize
## see limma user guide page 111 for example
## http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
setwd('dataExternal/gse29746/raw/')
untar('GSE29746_RAW.tar')
csFiles = list.files(pattern='*.gz', full.names = F)
oData = read.maimages(csFiles, source='agilent', green.only = T)
setwd(gcsWD)
# normalize the data
x.agilent = backgroundCorrect(oData, method='normexp')
x.agilent = normalizeBetweenArrays(x.agilent, method='quantile')
## create a new expression set object
# get matrix from normalised object
m = x.agilent@.Data[[1]]
# rename the probe ids
rownames(m) = x.agilent@.Data[[3]]$ProbeName
# match these row names to the GPL annotation
dfAn = fData(gse)
# sanity check
identical(as.character(dfAn$NAME), rownames(m))
# keep only the real genes
i = !is.na(dfAn$GENE)
table(i)
dfAn = dfAn[i,]; 
dfAn = droplevels.data.frame(dfAn)

m = m[i,]
# sanity check
identical(as.character(dfAn$NAME), rownames(m))
# remove duplicate probes/genes
i = !duplicated(rownames(m))
table(i)
m = m[i,]
dfAn = dfAn[i,]
dfAn = droplevels.data.frame(dfAn)
# sanity check
identical(as.character(dfAn$NAME), rownames(m))
# create expression set object
oExp = ExpressionSet(m)
fData(oExp) = dfAn

# get the samples from the expression set object
dfSamples = pData(gse)

# col names of the expression set object 
cn = colnames(exprs(oExp))
# remove the last .txt from the names
cn = gsub('^(GSM\\d+).txt', replacement = '\\1', cn, perl = T)

# order the sample names according to cel files
table(cn %in% rownames(dfSamples))
i = match(cn, rownames(dfSamples))
dfSamples = dfSamples[i,]
## sanity check
identical(rownames(dfSamples), cn)

## complete the expression set object by adding this sample information
pData(oExp) = dfSamples
colnames(exprs(oExp)) = cn
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

dbGetQuery(db, paste('describe MetaFile;'))
cFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

## create the entry for samples
cSampleCol
g_pid = 9
g_did = 19
title = paste(as.character(dfSamples$title), rownames(dfSamples))
description = paste(as.character(dfSamples$treatment_protocol_ch1))
group1 = rep('RA', times=nrow(dfSamples))
i = grep('HSF', as.character(dfSamples$characteristics_ch1.1))
group1[i] = 'Healthy'
i = grep('OASF', as.character(dfSamples$characteristics_ch1.1))
group1[i] = 'OA'

group2 = rep('female', times=nrow(dfSamples))
i = grep('Male', as.character(dfSamples$characteristics_ch1.2))
group2[i] = 'male'
group3 = gsub('\\D', '', as.character(dfSamples$characteristics_ch1.3))

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=title, description=description, group1=group1, 
                       group2=group2, group3=group3)
# write this data to the database
rownames(dfSamples) = NULL

# # write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)

## check annotation data
str(fData(oExp))
## store the grouping factors in affy object
oExp$fCondition = dfSamples$group1
oExp$fGender = dfSamples$group2
oExp$fAge = as.numeric(as.character(dfSamples$group3))

## check normalisation 
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

m = exprs(oExp)
range(m)
dim(m)
#m = m + abs(min(m)) + 0.1
oDiag.1 = CDiagnosticPlots(m, 'Raw Normalised format')

fBatch = oExp$fGender

## check normalisation
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)

plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)

plot.heatmap(oDiag.1)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

## this should give an error as scaling can't be done
## if all the vector 0 for PCA
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)

plot.PCA(oDiag.1.2, fBatch, legend.pos = 'topright')

plot.dendogram(oDiag.1.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## there seems to be a batch effect so use the batch as a covariate as well in model
h = hclust(dist(t(m)))
c = cutree(h, k = 3)
fBatch = factor(c)
oExp$fBatch = fBatch
### save the expression object and create metafile entry
getwd()
n = make.names(paste('GSE29746 Transcriptional analysis of normal and pathological synovial fibroblasts.rds'))
n2 = paste0('~/Data/MetaData/', n)
save(oExp, file=n2)

## access the mysql database to create entry for this data object
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/', comment='GSE29746 Transcriptional analysis of normal and pathological synovial fibroblasts. Raw data normalised using limma with green channel only and there is a batch effect which will be used as a covariate in the model.')
#dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

