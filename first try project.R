setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA")

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

# exprs functions returns the intensity values for each sample (column)
exprs(BreastCancer)
limma::plotDensities(exprs(BreastCancer)) # check this further 

## Reads in all .cel files and takes phenoData from the ExpressionFeatureset we loaded using https_ArrayExpress
BreastCancer <- ReadAffy(phenoData=pData(BreastCancer)) 

#pData returns the phenotypic data we loaded using ReadAffy
pData(BreastCancer)

## Quality Control on raw data
####################
# this is used to perform an explorative quality evaluation of hte orginal data and the log transformed data
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

# now can perform SAM or LIMMA to get the differentially expressed genes
#SAM
head(pData(BreastCancerRMA))
annot <- factor(pData(BreastCancerRMA)[,7]) # normal breast tissue and breast tumor tissue
length(annot) # 86

## Differential expression by SAM
annotb <- as.double(annot==annot[4]) # we want the breast tumor tissue to be one and the control to be zero 
sam.out_RMA <- sam(exprs(BreastCancerRMA),annotb)
summary(sam.out_RMA,3.1) # depends on you stringent you want to be 
# 559 identified genes with 0.03 falsely called genes or FDR of  2.95e-05
summary(sam.out_RMA,3.8)
# 273 identified genes with 0 FDR and 0 falsely called genes 
# the problem is we don't account for the persons 
pData(BreastCancerRMA)$Patients<- pData(BreastCancerRMA)$Hybridization.Name
pData(BreastCancerRMA)$Patients <- gsub('Normal', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub('Cancer', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub('T', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub('N', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients <- gsub(' ', '',pData(BreastCancerRMA)$Patients)
pData(BreastCancerRMA)$Patients # this is column with patients every patient has to occur 2 times in this column 
sum(grepl('BC0155', pData(BreastCancerRMA)$Patients))
sum(grepl('BC0117', pData(BreastCancerRMA)$Patients))
# zal wel kloppen

ID <- factor(pData(BreastCancerRMA)$Patients)
# vragen of we dit in ons model moeten opnemen wnt we missen dan wel veel vrijheidsgraden


## Differential expression by LIMMA
# Method as stated in limma package (no intercept, easy for simple model designs)
design <- model.matrix(~0+annot)
colnames(design)<-c("Cancer_breast_tissue","normal_breast_tissue")

fit <- lmFit(BreastCancerRMA,design)
cont.matrix <- makeContrasts(NvsS=Cancer_breast_tissue-normal_breast_tissue,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2)
volcanoplot(fit2)
limma::plotMA(fit2)
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(BreastCancerRMA)))
head(LIMMAout)
# To understand how this limma workflow works have a look at the bioconductor manual, the workflow with an 
# example is provided there
# Alternative (also applicable for more complex models)
design <- model.matrix(~annot)
fit <- lmFit(BreastCancerRMA,design)
fit2 <- eBayes(fit2)
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(MouseRMA)))
head(LIMMAout) 

# dit zijn dezelfde top genes als bij SAM!

## Check intensity values for top results
exprs(BreastCancerRMA)[rownames(exprs(BreastCancerRMA))%in%rownames(head(LIMMAout)),]
rowMeans(exprs(BreastCancerRMA)[rownames(exprs(BreastCancerRMA))%in%rownames(head(LIMMAout)),c(1,6,7)])
rowMeans(exprs(BreastCancerRMA)[rownames(exprs(BreastCancerRMA))%in%rownames(head(LIMMAout)),2:5])
# note that the probesets are not necessarly in the same order as for the limma output

## Optional: Differential expression analysis with MAS preprocessed data
####################
# we were working on the probe level not at level of genes per se
# most genes are represented by a single probeset on a microarray, for some genes multiple probesets are present. 

## Load annotation and sort alphabetically on probe name
annotation_BC <- read.table("A-AFFY-33.adf.txt",header=T,sep="\t",skip=17,fill=T)
print(head(annotation_BC))
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
sum(table(l*)) # no idea what this is for? 
LIMMAout_annot[LIMMAout_annot$gene=="ENSG00000119888",]
