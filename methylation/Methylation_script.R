setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/methylation/GSE101443_RAW")
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101443
#laptop 
setwd("~/documenten/Documenten/Ugent/applied high-throughput analysis/project_AHTA/methylation/GSE101443_RAW")

## Load packages
library('lumi')
library('wateRmelon')
library('ChAMPdata')


methyldata <- readEPIC(getwd())
ID <- c('GSM2703232','GSM2703233', 'GSM2703234', 'GSM2703235', 'GSM2703236','GSM2703237', 'GSM2703238', 'GSM2703239')
condition <- c('Tumour_A','Normal_A','Tumour_B','Normal_B','Tumour_C','Normal_C','Tumour_D','Normal_D')
Patient <- c('A','A','B','B','C','C','D','D')
annotation <- data.frame(ID,condition,Patient)

## Have a look at the data and annotation
print(methyldata)
print(dim(methyldata))
print(annotation)
print(sum(is.na(exprs(methyldata)))) # 1918
print(head(betas(methyldata))) # The “betas” function will retreive the beta values ( = methylation percentages) and the “exprs” function will retreive the M-values.
print(head(exprs(methyldata)))

## Change sampleNames to something more comprehensible
sampleNames(methyldata) <- annotation[,2]
sampleNames(methyldata)

## Remove NA values
methyldata <- methyldata[!(rowSums(is.na(exprs(methyldata)))>=1),]
methyldata

## Remove probes for which calling p-value is insufficient
methyldata.pf<-pfilter(methyldata) # removes the probes which we are not sure that are called correctly
#1408 probes were removed 

## Comparison of average methylation between control and tumor samples
boxplot(betas(methyldata),las=2)
meth_mean_tumour <- rep(0,8)
meth_mean_control <- rep(0,8)
for (i in 1:ncol(methyldata)){
  if((i %% 2) == 0) { # even 
    meth_mean_control[i] <- mean(betas(methyldata)[,i])
  } else {
    meth_mean_tumour[i] <- mean(betas(methyldata)[,i])
  }
}
meth_mean_tumour <- meth_mean_tumour[c(1,3,5,7)]
meth_mean_control <- meth_mean_control[c(2,4,6,8)]

t_test_res <- t.test(meth_mean_control,meth_mean_tumour,var.equal=F)
t_test_res

dat_boxplot <- data.frame(betas = c(meth_mean_tumour,meth_mean_control),
                          group = c('Tumour','Tumour','Tumour','Tumour','Normal','Normal','Normal','Normal'))
boxplot(betas~group,dat_boxplot,las=2)

## Normalization & QC
#######################

## Perform normalization including dye color adjustment
methyldata.dasen.pf <- dasen(methyldata.pf) 
#we will use the “dasen” function to adjust for color bias and normalize our data.

## Make methylumi objects to check density and color bias adjustment
#transform both the non-normalized counts as the normalized counts to a MethyLumiM object
methyldataM <- as(methyldata.pf, 'MethyLumiM') 
methyldataN <- as(methyldata.dasen.pf, 'MethyLumiM')

## Make QC plot
par(mfrow=c(2,1))
plotColorBias1D(methyldataM,channel="both",main="before",xlab="Beta-value")
plotColorBias1D(methyldataN,channel="both",main="after",xlab="Beta-value")
density(methyldataM,xlab="M-value",main="before")
density(methyldataN,xlab="M-value",main="after")

## Differential methylation analysis: limma
############################
par(mfrow=c(1,1))
## Build design and contrasts
condition <- factor(as.character(c('Tumour','Normal','Tumour','Normal','Tumour','Normal','Tumour','Normal')))
patient <- factor(annotation$Patient)

design2 <- model.matrix(~0+condition+patient)
colnames(design2)[1:2] <- c("Control","Tumor")
design2
cont.matrix2 <- makeContrasts(TumourvsControl=Tumor-Control,levels=design2)

## Limma
fit_2 <- lmFit(methyldataN,design2) # normalized data
fit_2 <- contrasts.fit(fit_2,cont.matrix2)
fit_2 <- eBayes(fit_2)
volcanoplot(fit2)
limma::plotMA(fit2, main= 'MA-plot')
LIMMAout_2 <- topTable(fit_2,adjust="BH",number=nrow(exprs(methyldata)))
LIMMAout_2

dim(LIMMAout_2[LIMMAout_2$adj.P.Val <= 0.0615913,]) 
dim(LIMMAout_2[abs(LIMMAout_2$logFC) > 2 & LIMMAout_2$adj.P.Val <= 0.0615913,]) # 358 significante genen gevonden 

dim(LIMMAout_2[LIMMAout_2$logFC > 2 & LIMMAout_2$adj.P.Val <= 0.0615913,]) # 251 significante genen gevonden 
dim(LIMMAout_2[LIMMAout_2$logFC < (-2) & LIMMAout_2$adj.P.Val <= 0.0615913,]) # 107 significante genen gevonden 
## Check M-values for top results
exprs(methyldataN)[rownames(methyldataN)%in%rownames(head(LIMMAout_2)),]
betas(methyldataN)[rownames(methyldataN)%in%rownames(head(LIMMAout_2)),]


## Functional annotation of limma results
############################

## Load annotation and sort alphabetically on probe name
#https://www.bioconductor.org/packages/devel/data/experiment/manuals/ChAMPdata/man/ChAMPdata.pdf
data("probe.features")

annotation_MA <- probe.features
print(head(annotation_MA))
annotation_MA <- annotation_MA[sort(rownames(annotation_MA),index.return=T)$ix,]

## Check if all probes are present in both sets
dim(LIMMAout_2)
dim(annotation_MA) # annotation has more rows than Limma output
sum(LIMMAout_2$Probe_ID%in%rownames(annotation_MA))
sum(rownames(annotation_MA)%in%LIMMAout_2$Probe_ID) 
# Also check the reverse so no duplicate rows are present in annotation

## Since more probes are present in the annotation file, remove unnecessary probes
annotation_MA <- annotation_MA[rownames(annotation_MA)%in%LIMMAout_2$Probe_ID,]
dim(annotation_MA)
## Sort LIMMA output alphabetically on probe name
LIMMAout_sorted_2 <- LIMMAout_2[sort(LIMMAout_2$Probe_ID,index.return=T)$ix,]

## Add gene names to LIMMA output
LIMMAout_sorted_2$Gene <- annotation_MA$gene
LIMMAout_sorted_2$Feature <- annotation_MA$feature
LIMMAout_sorted_2$Chrom <- annotation_MA$CHR
LIMMAout_sorted_2$Pos <- annotation_MA$MAPINFO
LIMMAout_sorted_2$Chrom <- as.character(LIMMAout_sorted$Chrom)
LIMMAout_sorted_2$Gene <- as.character(LIMMAout_sorted$Gene)
LIMMAout_sorted_2$Feature <- as.character(LIMMAout_sorted$Feature)

# The data type for these columns is altered to prevent issues further downstream
LIMMAout_annot_2 <- LIMMAout_sorted_2[sort(LIMMAout_sorted_2$P.Value,index.return=T)$ix,c(1,12,13,10,11,4,7,8)] 
LIMMAout_annot_2<- LIMMAout_annot_2[!LIMMAout_annot_2$Gene=="",]# filtering van de lege genes
LIMMAout_annot_2

dim(LIMMAout_annot_2[abs(LIMMAout_annot_2$logFC) > 2 & LIMMAout_annot_2$adj.P.Val <= 0.0615913,]) # we vinden 250 genen die significante differentiele methylatie hebben tumour vs control
dim(LIMMAout_annot_2[LIMMAout_annot_2$logFC > 2 & LIMMAout_annot_2$adj.P.Val <= 0.0615913,]) # 182
dim(LIMMAout_annot_2[LIMMAout_annot_2$logFC < (-2) & LIMMAout_annot_2$adj.P.Val <= 0.0615913,]) # 68
# een positive logFC wijst op hogere methylatie in tumour sample dan in de control
# een negative logFC wijst op een hogere methylatie in de control in vergelijking met de tumour.

significant_p_values <- LIMMAout_annot_2[abs(LIMMAout_annot_2$logFC) > 2 & LIMMAout_annot_2$adj.P.Val <= 0.0615913,]
foldchanges<- sort(significant_p_values$logFC,decreasing=TRUE)[0:10]

topten <- significant_p_values[significant_p_values$logFC%in%foldchanges,]
topten$Probe_ID
dim(significant_p_values)
# saving results
save(significant_p_values, file= "Methylation_significant.Rda")
load('Methylation_significant.Rda')
dim(significant_p_values)

pos_logfold<- significant_p_values[significant_p_values$logFC > 2 & significant_p_values$adj.P.Val <= 0.0615913,]
neg_logfold<- significant_p_values[significant_p_values$logFC < (-2) & significant_p_values$adj.P.Val <= 0.0615913,]


write.table(significant_p_values$Gene,file='methylation_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(pos_logfold$Gene,file='methylation_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$Gene,file='methylation_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
## Interpretation results
############################

## Select CpGs in genic regions

## Select CpGs in promoter regions
LIMMAout_annot_prom <- significant_p_values[grepl("TSS",significant_p_values$Feature) | (significant_p_values$Feature=="1stExon"),] #TSS = transcription start site
dim(LIMMAout_annot_prom) # 82 'significante genen' in de promoter regions 
head(LIMMAout_annot_prom)

write.table(LIMMAout_annot_prom$Gene,file='Promoter_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)



