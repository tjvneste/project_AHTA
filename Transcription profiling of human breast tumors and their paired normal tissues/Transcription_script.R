#Transcription profiling of human breast tumors and their paired normal tissues
#desktop
setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")
#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-15852/

"BACKGROUND
Microarray is widely used to monitor gene expression changes in breast cancer. The transcriptomic changes in breast cancer is commonly occured during the transition of normal cells to cancerous cells. 
This is the first study on gene expression profiling of multi ethnic of Malaysian breast cancer patients (Malays, Chinese and Indian). 
We aim to identify differentially expressed genes between tumors and normal tissues. We have identified a set of 33 significant differentially expressed genes in the tumor vs. normal group at p<0.001. 
We study the gene expression patterns of 43 breast tumors and their paired normal control by using Affymetrix genechip U133A. 
We have identified a set of 33 significant differentially expressed genes in the tumor vs. normal group at p<0.001. Experiment Overall Design: Total RNAs were extracted from breast cancer and normal tissues. 
Samples were processed and hybridized on the chip for 16 hours. 
At the end of the study, we obtained a total of 86 set of gene expression data, which were from 43 tumors and 43 normal tissues. The gene expression were then compared between the tumor and normal groups.
"

## Load packages
library(affy)
library(arrayQualityMetrics)
library(ArrayExpress)
library(limma)
library(siggenes)
library(biomaRt)

## Import Data
####################

## Load in the data 
## Download data to your working directory
getAE("E-GEOD-15852", type = 'raw')

# load in the expressionFeatureSet object
BreastCancer <- ArrayExpress("E-GEOD-15852")

# exprs functions returns the intensity values for each sample (column)
exprs(BreastCancer)

## Reads in all .cel files and takes phenoData from the ExpressionFeatureset we loaded using https_ArrayExpress
BreastCancer <- ReadAffy(phenoData=pData(BreastCancer)) 
BreastCancer <- ReadAffy() # bcs data won't load from arrayexpress

# Creating Annotation matrix 
pData(BreastCancer)
disease <- NULL
patients <- rep(1:43, each = 2)  


for (i in pData(BreastCancer)$sample){
  if((i %% 2) == 0) {
    disease[i] <- 'Cancer'}
  else{
    disease[i] <- 'Normal'}

}
disease
patients
annotation <- data.frame(annotation,patients)
dim(exprs(BreastCancer))


## Quality Control on raw data
####################
# this is used to perform an explorative quality evaluation of the orginal data and the log transformed data
## arrayQualityMetrics (open "index.html" file for a full overview of the output)
arrayQualityMetrics(BreastCancer,outdir="~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues/raw",force=T)
arrayQualityMetrics(BreastCancer,outdir="~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues/rawlog",force=T,do.logtransform=T)

# Preprocessing of the data 
# rma function is used to perform background correction as well as quantile normalization 
BreastCancerRMA<- affy::rma(BreastCancer,background=T) #  If TRUE, background correct using RMA background correction


## Quality Control on preprocessed data
## QC post preprocessing
arrayQualityMetrics(BreastCancerRMA,outdir="~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues/RMA",force=T)  			#RMA produces log-transformed data

head(exprs(BreastCancerRMA))
head(exprs(BreastCancer))
#annot <- factor(pData(BreastCancerRMA)[,7]) # normal breast tissue and breast tumor tissue
annotation
dim(annotation) # 86 2

#annotb <- as.double(annot==annot[4]) # we want the breast tumor tissue to be one and the control to be zero 
#annotb
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


annotation$patients <- factor(annotation$patients)
annotation$annotation <- factor(annotation$annotation)
annotation

## Differential expression by LIMMA
# why LIMMA? => Limma has more power than SAM
# Method as stated in limma package (no intercept, easy for simple model designs)
design <- model.matrix(~0+annotation+patients, data= annotation)
colSums(design) # checking if this is correct 
#patients toevoegen
colnames(design)[1:2]<-c("Cancer_tissue","normal_tissue")

fit <- lmFit(BreastCancerRMA,design)
cont.matrix <- makeContrasts(CancervsControl=Cancer_tissue-normal_tissue,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2)

par(mfrow=c(1,1))
volcanoplot(fit2)
limma::plotMA(fit2, main= 'MA-plot')

LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(BreastCancerRMA))) 
head(LIMMAout)

# Have a look at the results and search for other probesets for your DE genes
significant_pvalues<- LIMMAout[LIMMAout$adj.P.Val<0.05,]
dim(significant_pvalues) # 5333 probes statistical significant
significant_pvalues_1<- LIMMAout[LIMMAout$adj.P.Val<0.05 & abs(LIMMAout$logFC) >1,]
dim(significant_pvalues_1) # 40 probes statistical significant 
significant_pvalues_2<- LIMMAout[LIMMAout$adj.P.Val<0.05 & abs(LIMMAout$logFC) >1.5,]
dim(significant_pvalues_2) # 9 probes statistical significant 
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
significant_pvalues<- LIMMAout_annot[LIMMAout_annot$adj.P.Val<0.05& abs(LIMMAout$logFC) >1,]
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
biologisch_relevant <- significant_pvalues[abs(significant_pvalues$logFC)>1,]
