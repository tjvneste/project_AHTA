setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/miRNA profiling array results in serum RNA (healthy volunteers and breast tumor patients)")

library(affy)
library(arrayQualityMetrics)
library(ArrayExpress)
library(limma)
library(siggenes)
library(oligo)
## Import Data
####################

## Load in the data 
# If this dataset is an AffyBatch object you can proceed with the Quality control on the raw data. If this 
# is an ExpressionFeatureSet instead continue to get the right object. We need an Affybatch object to secure
# compatibility with other packages used.

## Download data to your working directory
getAE("E-MTAB-6652")

dat <- read.celfiles(list.celfiles())
dat

BreastCancer_miRNA <- ReadAffy(phenoData=pData(dat)) 

# exprs functions returns the intensity values for each sample (column)
exprs(dat)
dim(exprs(dat))
#pData returns the phenotypic data we loaded using ReadAffy
pData(dat)

## Quality Control on raw data
####################
# this is used to perform an explorative quality evaluation of hte orginal data and the log transformed data
## arrayQualityMetrics (open "index.html" file for a full overview of the output)
arrayQualityMetrics(dat,outdir="/Users/tristanvanneste/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/miRNA profiling array results in serum RNA (healthy volunteers and breast tumor patients)/raw",force=T)
arrayQualityMetrics(dat,outdir="/Users/tristanvanneste/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/miRNA profiling array results in serum RNA (healthy volunteers and breast tumor patients)/rawlog",force=T,do.logtransform=T)

# Preprocessing of the data 
# rma function is used to perform background correction as well as quantile normalization 
BreastCancer_miRNA <- rma(dat)
BreastCancer_miRNA
limma::plotDensities(exprs(BreastCancer_miRNA))
## Quality Control on preprocessed data
## QC post preprocessing
arrayQualityMetrics(BreastCancer_miRNA,outdir="/Users/tristanvanneste/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/miRNA profiling array results in serum RNA (healthy volunteers and breast tumor patients)/miRNA",force=T)  			#RMA produces log-transformed data

# now can perform SAM or LIMMA to get the differentially expressed genes
#SAM
pData(BreastCancer_miRNA)
exprs(BreastCancer_miRNA)
annot <- factor(pData(BreastCancer_miRNA)) # normal breast tissue and breast tumor tissue
length(annot) # 86

## Differential expression by SAM
annotb <- as.double(annot==annot[4]) # we want the breast tumor tissue to be one and the control to be zero 
sam.out_RMA <- sam(exprs(BreastCancer_miRNA),annotb)
summary(sam.out_RMA,3.1) # depends on you stringent you want to be 
# 559 identified genes with 0.03 falsely called genes or FDR of  2.95e-05
summary(sam.out_RMA,3.8)
# 273 identified genes with 0 FDR and 0 falsely called genes 
# the problem is we don't account for the persons 
pData(BreastCancer_miRNA)$Patients<- pData(BreastCancer_miRNA)$Hybridization.Name
pData(BreastCancer_miRNA)$Patients <- gsub('Normal', '',pData(BreastCancer_miRNA)$Patients)
pData(BreastCancer_miRNA)$Patients <- gsub('Cancer', '',pData(BreastCancer_miRNA)$Patients)
pData(BreastCancer_miRNA)$Patients <- gsub('T', '',pData(BreastCancer_miRNA)$Patients)
pData(BreastCancer_miRNA)$Patients <- gsub('N', '',pData(BreastCancer_miRNA)$Patients)
pData(BreastCancer_miRNA)$Patients <- gsub(' ', '',pData(BreastCancer_miRNA)$Patients)
pData(BreastCancer_miRNA)$Patients # this is column with patients every patient has to occur 2 times in this column 
sum(grepl('BC0155', pData(BreastCancer_miRNA)$Patients))
sum(grepl('BC0117', pData(BreastCancer_miRNA)$Patients))
# zal wel kloppen

ID <- factor(pData(BreastCancer_miRNA)$Patients)
# vragen of we dit in ons model moeten opnemen wnt we missen dan wel veel vrijheidsgraden


## Differential expression by LIMMA
# Method as stated in limma package (no intercept, easy for simple model designs)
design <- model.matrix(~0+annot)
colnames(design)<-c("Cancer_breast_tissue","normal_breast_tissue")

fit <- lmFit(BreastCancer_miRNA,design)
cont.matrix <- makeContrasts(NvsS=Cancer_breast_tissue-normal_breast_tissue,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2)
volcanoplot(fit2)
limma::plotMA(fit2)
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(BreastCancer_miRNA)))
head(LIMMAout)
# To understand how this limma workflow works have a look at the bioconductor manual, the workflow with an 
# example is provided there
# Alternative (also applicable for more complex models)
design <- model.matrix(~annot)
fit <- lmFit(BreastCancer_miRNA,design)
fit2 <- eBayes(fit2)
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(MouseRMA)))
head(LIMMAout) 

# dit zijn dezelfde top genes als bij SAM!

## Check intensity values for top results
exprs(BreastCancer_miRNA)[rownames(exprs(BreastCancer_miRNA))%in%rownames(head(LIMMAout)),]
rowMeans(exprs(BreastCancer_miRNA)[rownames(exprs(BreastCancer_miRNA))%in%rownames(head(LIMMAout)),c(1,6,7)])
rowMeans(exprs(BreastCancer_miRNA)[rownames(exprs(BreastCancer_miRNA))%in%rownames(head(LIMMAout)),2:5])
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
