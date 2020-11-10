setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/methylation")


library(ArrayExpress)

## Download data to your working directory
getAE("E-GEOD-66695")

methylation <- ArrayExpress("E-GEOD-66695")
methylation <- BreastCancer
