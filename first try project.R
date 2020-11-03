setwd("~/Documents/Bioinformatics/project AHATA")

BiocManager::install("affy",update=F)
BiocManager::install("arrayQualityMetrics",update=F)
BiocManager::install("ArrayExpress",update=F)
BiocManager::install("limma",update=F)
BiocManager::install("siggenes",update=F)

## Load packages
library(affy)
library(arrayQualityMetrics)
library(ArrayExpress)
library(limma)
library(siggenes)

## Import Data
####################

## Load in the data 
# If this dataset is an AffyBatch object you can proceed with the Quality control on the raw data. If this 
# is an ExpressionFeatureSet instead continue to get the right object. We need an Affybatch object to secure
# compatibility with other packages used.

## Download data to your working directory
getAE("E-GEOD-15852", type = 'raw')

# load in the expressionFeatureSet object
BreastCancer <- ArrayExpress("E-GEOD-15852")
BreastCancer

## Reads in all .cel files and takes phenoData from the ExpressionFeatureset we loaded using https_ArrayExpress
BreastCancer <- ReadAffy(phenoData=pData(BreastCancer)) 

## Have a look to the data you just loaded
head(exprs(BreastCancer))
pData(BreastCancer)

## Quality Control on raw data
####################

## arrayQualityMetrics (open "index.html" file for a full overview of the output)
arrayQualityMetrics(BreastCancer,outdir="/Users/tristanvanneste/Documents/Bioinformatics/project Tristan/raw",force=T)
arrayQualityMetrics(BreastCancer,outdir="/Users/tristanvanneste/Documents/Bioinformatics/project Tristan/rawlog",force=T,do.logtransform=T)


