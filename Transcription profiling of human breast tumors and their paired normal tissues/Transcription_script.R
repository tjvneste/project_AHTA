#Transcription profiling of human breast tumors and their paired normal tissues

setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")
#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-15852/
BiocManager::install("affy",update=F)
## Load packages
library(affy)
library(arrayQualityMetrics)
library(ArrayExpress)
library(limma)
library(siggenes)
library(ggplot2)
library(biomaRt)

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

# exprs functions returns the intensity values for each sample (column)
exprs(BreastCancer)

## Reads in all .cel files and takes phenoData from the ExpressionFeatureset we loaded using https_ArrayExpress
BreastCancer <- ReadAffy(phenoData=pData(BreastCancer)) 

#pData returns the phenotypic data we loaded using ReadAffy
pData(BreastCancer)

## Quality Control on raw data
####################
# this is used to perform an explorative quality evaluation of the orginal data and the log transformed data
## arrayQualityMetrics (open "index.html" file for a full overview of the output)
arrayQualityMetrics(BreastCancer,outdir="/Users/tristanvanneste/Documents/Bioinformatics/project Tristan/raw",force=T)
arrayQualityMetrics(BreastCancer,outdir="/Users/tristanvanneste/Documents/Bioinformatics/project Tristan/rawlog",force=T,do.logtransform=T)

# Preprocessing of the data 
# rma function is used to perform background correction as well as quantile normalization 
BreastCancerRMA<- affy::rma(BreastCancer,background=T)
limma::plotDensities(exprs(BreastCancerRMA))
## Quality Control on preprocessed data
## QC post preprocessing
arrayQualityMetrics(BreastCancerRMA,outdir="/Users/tristanvanneste/Documents/Bioinformatics/project Tristan/RMA",force=T)  			#RMA produces log-transformed data

head(pData(BreastCancerRMA))
annot <- factor(pData(BreastCancerRMA)[,7]) # normal breast tissue and breast tumor tissue
annot
length(annot) # 86

## Differential expression by SAM
annotb <- as.double(annot==annot[4]) # we want the breast tumor tissue to be one and the control to be zero 
annotb
# the problem is we don't account for the persons but can we? bcs we have 86 samples? YES
pData(BreastCancerRMA)$Patients<- pData(BreastCancerRMA)$Hybridization.Name
pData(BreastCancerRMA)$Patients <- gsub('Normal', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub('Cancer', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub('T', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub('N', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub(' ', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients # this is column with patients every patient has to occur 2 times in this column 
sum(grepl('BC0155', pData(BreastCancerRMA)$Patients))
sum(grepl('BC0117', pData(BreastCancerRMA)$Patients))
ID <- factor(pData(BreastCancerRMA)$Patients)
length(levels(ID))

## Differential expression by LIMMA
# why LIMMA? => Limma has more power than SAM
# Method as stated in limma package (no intercept, easy for simple model designs)
design <- model.matrix(~0+annot+ID)
colSums(design) # checking if this is correct 
#patients toevoegen
colnames(design)[1:2]<-c("Cancer_breast_tissue","normal_breast_tissue")

fit <- lmFit(BreastCancerRMA,design)
cont.matrix <- makeContrasts(NvsS=Cancer_breast_tissue-normal_breast_tissue,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2)
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(BreastCancerRMA))) 
head(LIMMAout)

# Have a look at the results and search for other probesets for your DE genes
head(LIMMAout)
significant_pvalues_1<- LIMMAout[LIMMAout$adj.P.Val<0.05,]
dim(significant_pvalues_1)
head(significant_pvalues_1)
tail(significant_pvalues_1)

## Optional: Differential expression analysis with MAS preprocessed data
####################
# we were working on the probe level not at level of genes per se
# most genes are represented by a single probeset on a microarray, for some genes multiple probesets are present. 

## Load annotation and sort alphabetically on probe name
annotation_BC <- read.table("A-AFFY-33.adf.txt",header=T,sep="\t",skip=17,fill=T)
print(tail(annotation_BC))
annotation_BC[100:105,]
annotation_BC <- annotation_BC[sort(annotation_BC$Composite.Element.Name,index.return=T)$ix,]

## Check if all probes are present in both sets
dim(annotation_BC)
dim(LIMMAout)

## Double check => "Assumption is the mother of all fuck up's ;)"
sum(annotation_BC$Composite.Element.Name==sort(rownames(LIMMAout)))

## Sort LIMMA output alphabetically on probe name
LIMMAout_sorted <- LIMMAout[sort(rownames(LIMMAout),index.return=T)$ix,]
dim(LIMMAout_sorted)
## Add gene names to LIMMA output
LIMMAout_sorted$gene <- annotation_BC$Composite.Element.Database.Entry.ensembl.
LIMMAout_annot <- LIMMAout_sorted[sort(LIMMAout_sorted$adj.P.Val,index.return=T)$ix,]



# Have a look at the results and search for other probesets for your DE genes
head(LIMMAout_annot)
significant_pvalues<- LIMMAout_annot[LIMMAout_annot$adj.P.Val<0.05,]
dim(significant_pvalues)
head(significant_pvalues)
tail(significant_pvalues)

affyids <- rownames(significant_pvalues)

output_sign <-getBM(attributes = c('affy_hg_u133_plus_2', 'entrezgene_id','hgnc_symbol','chromosome_name', 'start_position', 'end_position'), #  is a vector of attributes that one wants to retrieve (= the output of the query).
      filters = 'affy_hg_u133_plus_2', #  is a vector of filters that one wil use as input to the query.
      values = affyids, # a vector of values for the filters
      mart = ensembl)


significant_pvalues[,"hgnc_symbol"] <-NA
significant_pvalues[,"chromosome_name"] <-NA
significant_pvalues[,"start_position"] <-NA
significant_pvalues[,"end_position"] <-NA
# filtering of the zero entries 
output_sign <- output_sign[!output_sign$hgnc_symbol=="",]
output_sign <- output_sign[!output_sign$affy_hg_u133_plus_2=="",]
count=1
significant_pvalues <- significant_pvalues[rownames(significant_pvalues)%in%output_sign$affy_hg_u133_plus_2,]

for (i in rownames(significant_pvalues)){
  significant_pvalues$hgnc_symbol[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,3]
  significant_pvalues$chromosome_name[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,4]
  significant_pvalues$start_position[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,5]
  significant_pvalues$end_position[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,6]
  count=count+1
}

save(significant_pvalues,file="Significant_output_annotation_transcription.Rda") # object noemt significant_pvalues
load("Significant_output_annotation_transcription.Rda")
