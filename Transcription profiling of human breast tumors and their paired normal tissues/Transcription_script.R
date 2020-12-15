#set working directory 
#https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-15852/
setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")

## Load packages
library(affy)
library(arrayQualityMetrics)
library(ArrayExpress)
library(limma)
library(biomaRt)

## Import Data
####################

## Load in the data 
## Download data to your working directory
getAE("E-GEOD-15852", type = 'raw')

# load in the expressionFeatureSet object
BreastCancer <- ArrayExpress("E-GEOD-15852")

## Reads in all .cel files and takes phenoData from the ExpressionFeatureset we loaded using https_ArrayExpress
BreastCancer <- ReadAffy(phenoData=pData(BreastCancer)) 
BreastCancer <- ReadAffy() # If ArrayExpress does not work => use this

dim(pData(BreastCancer))
exprs(BreastCancer)
pData(BreastCancer)

# Creating Annotation matrix 

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
annotation <- data.frame(disease,patients)
annotation
dim(exprs(BreastCancer))


## Quality Control on raw data
####################
# raw data
arrayQualityMetrics(BreastCancer,outdir="~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/
                    Transcription profiling of human breast tumors and their paired normal tissues/raw",force=T)
# logtransformed data
arrayQualityMetrics(BreastCancer,outdir="~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/
                    Transcription profiling of human breast tumors and their paired normal tissues/rawlog",force=T,do.logtransform=T)

# Preprocessing of the data 
BreastCancerRMA<- affy::rma(BreastCancer,background=T)


## Quality Control on preprocessed data
## QC post preprocessing
arrayQualityMetrics(BreastCancerRMA,outdir="~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/
                    Transcription profiling of human breast tumors and their paired normal tissues/RMA",force=T)  			#RMA produces log-transformed data

head(exprs(BreastCancerRMA))
head(exprs(BreastCancer))
#annot <- factor(pData(BreastCancerRMA)[,7]) # normal breast tissue and breast tumor tissue
annotation
dim(annotation) # 86 2

# this piece is only necessary if you got the object of ArrayExpress otherwise skip this
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

# begin here again 
annotation$patients <- factor(annotation$patients)
annotation$disease <- factor(annotation$disease)
annotation

## Differential expression by LIMMA
design <- model.matrix(~0+disease+patients, data= annotation)
colSums(design)
colnames(design)[1:2]<-c("Cancer_tissue","normal_tissue")

fit <- lmFit(BreastCancerRMA,design)
cont.matrix <- makeContrasts(CancervsControl=Cancer_tissue-normal_tissue,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2)

LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(BreastCancerRMA))) 
head(LIMMAout)

# MA plot
threshold.sign <- LIMMAout[LIMMAout$adj.P.Val<0.05,]
dim(threshold.sign)
with(LIMMAout, plot(AveExpr, logFC, pch=20,main="MA plot"))
with(subset(threshold.sign),points(AveExpr, logFC,pch=20,col="red"))
# histogram p-values
hist(fit2$p.value, main= 'distributions of the p-values',xlab='p-values')
# Volcano plot
threshold.sign <- LIMMAout[LIMMAout$adj.P.Val<0.05,]
dim(threshold.sign)
with(LIMMAout, plot(logFC, -log10(P.Value), pch=20,main="Volcano plot"))
with(subset(threshold.sign),points(logFC, -log10(P.Value),pch=20,col="red"))

     
# Have a look at the results 
significant_pvalues<- LIMMAout[LIMMAout$adj.P.Val<0.05,]
dim(significant_pvalues) # 5333 probes statistical significant
significant_pvalues_1<- LIMMAout[LIMMAout$adj.P.Val<0.05 & abs(LIMMAout$logFC) >1,]
dim(significant_pvalues_1) # 40 probes statistical significant 
significant_pvalues_1

significant_pvalues_2<- LIMMAout[LIMMAout$adj.P.Val<0.05 & abs(LIMMAout$logFC) >1.5,]
dim(significant_pvalues_2) # 9 probes statistical significant 
head(significant_pvalues_1)

## Load annotation and sort alphabetically on probe name
annotation_BC <- read.table("A-AFFY-33.adf.txt",header=T,sep="\t",skip=17,fill=T)
annotation_BC <- annotation_BC[sort(annotation_BC$Composite.Element.Name,index.return=T)$ix,]

## Check if all probes are present in both sets
dim(annotation_BC)
dim(LIMMAout)

## Double check => "Assumption is the mother of all fuck up's ;)"
sum(annotation_BC$Composite.Element.Name==sort(rownames(LIMMAout)))

## Sort LIMMA output alphabetically on probe name
LIMMAout_sorted <- LIMMAout[sort(rownames(LIMMAout),index.return=T)$ix,]
## Add gene names to LIMMA output
LIMMAout_sorted$gene <- annotation_BC$Composite.Element.Database.Entry.ensembl.
LIMMAout_annot <- LIMMAout_sorted[sort(LIMMAout_sorted$adj.P.Val,index.return=T)$ix,]



# Have a look at the results and search for other probesets for your DE genes
head(LIMMAout_annot)
dim(LIMMAout_annot)
significant_pvalues<- LIMMAout_annot[LIMMAout_annot$adj.P.Val<0.05& abs(LIMMAout_annot$logFC) >1,]
adjusted_pvalues<- LIMMAout_annot[LIMMAout_annot$adj.P.Val<0.05,]
dim(significant_pvalues) # 40
dim(adjusted_pvalues) # 5333 7

affyids <- rownames(adjusted_pvalues)

ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
output_sign <-getBM(attributes = c('affy_hg_u133_plus_2', 'entrezgene_id','hgnc_symbol','chromosome_name', 'start_position', 'end_position'),
                    filters = 'affy_hg_u133_plus_2', #  is a vector of filters that one will use as input to the query.
                    values = affyids, # a vector of values for the filters
                    mart = ensembl)

adjusted_pvalues[,"hgnc_symbol"] <-NA
adjusted_pvalues[,"chromosome_name"] <-NA
adjusted_pvalues[,"start_position"] <-NA
adjusted_pvalues[,"end_position"] <-NA
# filtering of the zero entries 

output_sign <- output_sign[!output_sign$hgnc_symbol=="",]
output_sign <- output_sign[!output_sign$affy_hg_u133_plus_2=="",]
count=1
adjusted_pvalues <- adjusted_pvalues[rownames(adjusted_pvalues)%in%output_sign$affy_hg_u133_plus_2,]

for (i in rownames(adjusted_pvalues)){
  adjusted_pvalues$hgnc_symbol[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,3]
  adjusted_pvalues$chromosome_name[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,4]
  adjusted_pvalues$start_position[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,5]
  adjusted_pvalues$end_position[count] <- output_sign[output_sign$affy_hg_u133_plus_2==i,6]
  count=count+1
}

dim(adjusted_pvalues) # 5028 11
head(adjusted_pvalues)

save(adjusted_pvalues,file="adjusted_pvalues_annotation_transcription.Rda") # object = adjusted_pvalues



# for the abs(logFC) > 1
affyids <- rownames(significant_pvalues)

ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
output_sign <-getBM(attributes = c('affy_hg_u133_plus_2', 'entrezgene_id','hgnc_symbol','chromosome_name', 'start_position', 'end_position'), 
                    filters = 'affy_hg_u133_plus_2', #  is a vector of filters that one will use as input to the query.
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
significant_pvalues2 <- significant_pvalues
dim(significant_pvalues2) # 39 11
head(significant_pvalues2)

save(significant_pvalues2,file="Significant_output_annotation_transcription.Rda") # object = significant_pvalues2
sessionInfo()
