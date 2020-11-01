#!usr/bin/Rscript

# ArrayID = E-GEOD-45666 ()

setwd("C:/Users/Boris/Documents/School/2de master/Applied high throughput analysis/Project/ExpressionArray/E-GEOD-45666.raw.1/") # boris 
setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA") # Tristan 

BiocManager::install("GEOquery")
#BiocManager::install("ArrayExpress")
#BiocManager::install("arrayQualityMetrics")

library("ArrayExpress")
library("arrayQualityMetrics")
library("limma")
library("GEOquery")

#download the data to your working directory
gds <- getGEO("GSE45666")
gsm <- getGEO(filename="GSE45666_RAW")
GSE45666 <- getGEO('GSE45666',GSEMatrix=TRUE)
show(GSE45666)

getAE("E-GEOD-45666")
# load in the expressionFeatureSet object
BreastCancer_exprArray_AE <- ArrayExpress("E-GEOD-45666")
dim(pData(BreastCancer_exprArray_AE)) # 116 52
# Load in the Affybatch object
BreastCancer_exprArray <- ReadAffy(phenoData=pData(BreastCancer_exprArray_AE))

# limma package, zal wss beter lukken dan agilp
targets <- readTargets("Targets.txt")
head(targets)

targets$DiseaseStatus<- as.factor(targets$DiseaseStatus)
summary(targets$DiseaseStatus)

x <- read.maimages(targets$FileName, source="agilent", green.only=TRUE)
dim(x)
# Om annotatie te laden, de correcte array moet wel ingevuld worden
require(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'wikigene_description',
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name'))

write.table(
  annotLookup,
  paste0('Agilent-021827 Human miRNA Microarray G4470C', gsub("-", "_", as.character(Sys.Date())), '.tsv'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE)


