# Name: 07_diagnosticsDataSets.R
# Auth: umar niazi umar.niazi@kcl.ac.uk
# Date: 25/10/2017
# Desc: perform diagnostics for each downloaded dataset

source('header.R')

library(Biobase)
library(downloader)

# load data from database
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and metafile table
dbGetQuery(db, paste('describe MetaFile;'))
cMetaFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

q = paste('select MetaFile.*, Data.id as did, Data.idProject as pid from MetaFile, Data
          where Data.idProject = 9 and Data.id = MetaFile.idData;')

dfMetaFile = dbGetQuery(db, q)

## produce report on each data set
n = paste0(dfMetaFile$location, dfMetaFile$name)

lExp = lapply(n, f_LoadObject) 

cn = lapply(lExp, function(x){
  # get the column names with factors
  # these should be the last 3 columns
  c = length(colnames(pData(x)))
  return(c(c-2, c-1, c))
})

lBatch = lapply(lExp, function(x){
  # get the main factor
  # these should be one of the last 3 columns
  c = grep('fCondition', colnames(pData(x)))
  return(pData(x)[,c])
})

dbDisconnect(db)
## make normalisation plots
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

names(lExp) = c('GSE7669', 'GSE56409', 'GSE83147', 'GSE29746', 'GSE49604.GPL10558', 'GSE49604.GPL8432')

lOb = lapply(seq_along(lExp), function(dat){
  o = CDiagnosticPlots(exprs(lExp[[dat]]), csTitle = names(lExp)[[dat]])
  l = CDiagnosticPlotsGetParameters(o)
  l$PCA.jitter = F
  l$HC.jitter = F
  CDiagnosticPlotsSetParameters(o, l)
})

sapply(seq_along(lOb), function(x) boxplot.median.summary(lOb[[x]], lBatch[[x]]))
sapply(seq_along(lOb), function(x) plot.PCA(lOb[[x]], lBatch[[x]]))
sapply(seq_along(lOb), function(x) plot.dendogram(lOb[[x]], lBatch[[x]], labels_cex = 0.7))

