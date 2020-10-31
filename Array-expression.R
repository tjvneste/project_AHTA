#!usr/bin/Rscript

# ArrayID = E-GEOD-45666 ()

setwd("C:/Users/Boris/Documents/School/2de master/Applied high throughput analysis/Project/ExpressionArray/E-GEOD-45666.raw.1/")

#BiocManager::install("agilp")
library("ArrayExpress")
library("arrayQualityMetrics")
library("limma")

#getAE("E-GEOD-45666")
#BreastCancer_exprArray_AE <- ArrayExpress("E-GEOD-45666")

# limma package, zal wss beter lukken dan agilp
targets <- readTargets("Targets.txt")
x <- read.maimages(targets$FileName, source="agilent", green.only=TRUE)

# Om annotatie te laden, de correcte array moet wel ingevuld worden
require(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'Agilent-021827 Human miRNA Microarray G4470C',
    'wikigene_description',
    'ensembl_gene_id',
    'entrezgene',
    'gene_biotype',
    'external_gene_name'))
write.table(
  annotLookup,
  paste0('Agilent-021827 Human miRNA Microarray G4470C', gsub("-", "_", as.character(Sys.Date())), '.tsv'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE)


