# Name: 01_downloadData_GSE7669.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 19/10/2017
# Desc: Download the relevant datasets, normalise and create database entries

source('header.R')

library(GEOquery)
library(downloader)

## open the soft format and raw data
gse =  getGEO(filename = 'dataExternal/gse7669/GSE7669_series_matrix.txt.gz')

library(affy)
setwd('dataExternal/gse7669/')
untar('GSE7669_RAW.tar')
oData = ReadAffy()
setwd(gcsWD)
# normalize the data
x.affy = rma(oData)

# get the samples from the expression set object
dfSamples = pData(gse)

# col names of the affy expression matrix
cn = colnames(exprs(x.affy))
# remove the last .CEL.gz from the names
cn = gsub('^(GSM\\d+).CEL\\.gz', replacement = '\\1', cn, perl = T)

# order the sample names according to cel files
table(cn %in% rownames(dfSamples))
i = match(cn, rownames(dfSamples))
dfSamples = dfSamples[i,]

## complete the expression set object
pData(x.affy) = dfSamples
colnames(exprs(x.affy)) = cn
# sanity check
t1 = colnames(x.affy)
t2 = rownames(pData(x.affy))
table(t1 %in% t2)
all(t1 == t2)
## create grouping factors
cvSource = as.character(pData(x.affy)[,'source_name_ch1'])
table(cvSource)
## create factors for cell type
fCell.type = rep(NA, length=length(cvSource))
i = grep('fibroblast', cvSource, ignore.case = T)
fCell.type[i] = 'fibroblast'
fCell.type[-i] = 'keratinocyte'
fCell.type = factor(fCell.type)
## keloid vs normal individual
fDisorder = rep(NA, length=length(cvSource))
i = grep('keloid', cvSource, ignore.case = T)
fDisorder[i] = 'keloid'
i = grep('without keloid', cvSource, ignore.case = T)
fDisorder[i] = 'healthy'
i = grep('non-lesion', cvSource, ignore.case = T)
fDisorder[i] = 'keloid.non-lesion'
# sanity check
data.frame(fDisorder, cvSource)
fDisorder = factor(fDisorder, levels = c('keloid', 'keloid.non-lesion', 'healthy'))
table(fDisorder); levels(fDisorder)
x.affy$fDisorder = fDisorder
x.affy$fCell.type = fCell.type










url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(d, 'method 1')

fBatch = as.character(oExp$title)
fBatch = factor(gsub('fibroblasts_(\\w+)_patie.+', '\\1', fBatch))

## compare the 2 methods using various plots
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)

plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## how about the extreme values
oDiag.1 = calculateExtremeValues(oDiag.1)
oDiag.2 = calculateExtremeValues(oDiag.2)
m1 = mGetExtremeValues(oDiag.1)
m2 = mGetExtremeValues(oDiag.2)

## samples with most extreme values
apply(m1, 2, function(x) sum(x > 0))
apply(m2, 2, function(x) sum(x > 0))

## variables that are contributing to this
v1 = apply(m1, 1, function(x) sum(x > 0))
v2 = apply(m2, 1, function(x) sum(x > 0))

which(v1 > 0)
which(v2 > 0)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

## this should give an error as scaling can't be done
## if all the vector 0 for PCA
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)
## reset flag for jittering
l$PCA.jitter = T
## should work this time
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)
plot.PCA(oDiag.1, fBatch, legend.pos = 'topright')
plot.PCA(oDiag.1.2, fBatch)
plot.dendogram(oDiag.1, fBatch)
plot.dendogram(oDiag.1.2, fBatch)


##############################################################################
