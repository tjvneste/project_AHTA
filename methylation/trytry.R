setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/methylation/GSE101443_RAW")

## Load packages
library('lumi')
library('wateRmelon')
library('ChAMPdata')


methyldata <- readEPIC(getwd())
ID <- c('GSM2703232','GSM2703233', 'GSM2703234', 'GSM2703235', 'GSM2703236','GSM2703237', 'GSM2703238', 'GSM2703239')
condition <- c('Tumour_A','Normal_A','Tumour_B','Normal_B','Tumour_C','Normal_C','Tumour_D','Normal_D')
annotation <- data.frame(ID,condition)

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
methyldata.pf<-pfilter(methyldata) 

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

dat_boxplot <- data.frame(betas = c(meth_mean_control,meth_mean_control),
                          group = c('Tumour','Normal','Tumour','Normal','Tumour','Normal','Tumour','Normal'))
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
par(mfrow=c(2,2))
plotColorBias1D(methyldataM,channel="both",main="before")
plotColorBias1D(methyldataN,channel="both",main="after")
density(methyldataM,xlab="M-value",main="before")
density(methyldataN,xlab="M-value",main="after")

## Differential methylation analysis: limma
############################

## Build design and contrasts
des <- factor(as.character(c('Tumour','Normal','Tumour','Normal','Tumour','Normal','Tumour','Normal')))
design <- model.matrix(~0+des)
colnames(design) <- c("Normal","Tumor")
cont.matrix <- makeContrasts(NvsS=Tumor-Normal,levels=design)

## Limma
fit <- lmFit(methyldataN,design) # normalized data
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
volcanoplot(fit2)
limma::plotMA(fit2)
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(methyldata)))
head(LIMMAout)

## Check M-values for top results
exprs(methyldataN)[rownames(methyldataN)%in%rownames(head(LIMMAout)),]
# note that the probesets are not necessarly in the same order as for the limma output

## Functional annotation of limma results
############################

## Load annotation and sort alphabetically on probe name
#https://www.bioconductor.org/packages/devel/data/experiment/manuals/ChAMPdata/man/ChAMPdata.pdf
#data("probe.features.epic")
data("probe.features") # this one?  

annotation_MA <- probe.features
print(head(annotation_MA))
annotation_MA <- annotation_MA[sort(rownames(annotation_MA),index.return=T)$ix,]

## Check if all probes are present in both sets
dim(LIMMAout)
dim(annotation_MA)
sum(LIMMAout$Probe_ID%in%rownames(annotation_MA))
sum(rownames(annotation_MA)%in%LIMMAout$Probe_ID) 
# Also check the reverse so no duplicate rows are present in annotation

## Since more probes are present in the annotation file, remove unnecessary probes
annotation_MA <- annotation_MA[rownames(annotation_MA)%in%LIMMAout$Probe_ID,]

## Sort LIMMA output alphabetically on probe name
LIMMAout_sorted <- LIMMAout[sort(LIMMAout$Probe_ID,index.return=T)$ix,]

## Add gene names to LIMMA output
LIMMAout_sorted$Gene <- annotation_MA$gene
LIMMAout_sorted$Feature <- annotation_MA$feature
LIMMAout_sorted$Chrom <- annotation_MA$CHR
LIMMAout_sorted$Pos <- annotation_MA$MAPINFO
LIMMAout_sorted$Chrom <- as.character(LIMMAout_sorted$Chrom)
LIMMAout_sorted$Gene <- as.character(LIMMAout_sorted$Gene)
LIMMAout_sorted$Feature <- as.character(LIMMAout_sorted$Feature)
# The data type for these columns is altered to prevent issues further downstream
LIMMAout_annot <- LIMMAout_sorted[sort(LIMMAout_sorted$P.Value,index.return=T)$ix,c(1,12,13,10,11,4,7,8,5)] 
# Sort on p-values to prevent errors in sorting due to equal FDR values


## Interpretation results
############################

## Select CpGs in genic regions
sum(LIMMAout_annot$adj.P.Val<0.05)
sum(LIMMAout_annot$adj.P.Val[LIMMAout_annot$Gene!=""]<0.05)

LIMMAout_annot_gene <- LIMMAout_annot[LIMMAout_annot$Gene!="",]

## Check genic results 
head(LIMMAout_annot_gene)
LIMMAout_annot_gene[0:50,]
LIMMAout_annot_gene[50:100,]
topgenes_genic <- unique(LIMMAout_annot_gene$Gene[1:10])
for (i in 1:length(topgenes_genic)){
  LIMMAout_subset <- LIMMAout_annot_gene[(LIMMAout_annot_gene$Gene==topgenes_genic[i]) & 
                                           (LIMMAout_annot_gene$adj.P.Val<0.05) & (abs(LIMMAout_annot_gene$logFC)>2),]
  print(LIMMAout_subset[sort(LIMMAout_subset$Pos,index.return=T)$ix,])
}

## Select CpGs in promoter regions
LIMMAout_annot_prom <- LIMMAout_annot_gene[grepl("TSS",LIMMAout_annot_gene$Feature) | (LIMMAout_annot_gene$Feature=="1stExon"),]
head(LIMMAout_annot_prom)

## Look for multiple CpG in promoter regions undergoing similar methylation differences
topgenes_prom <- unique(LIMMAout_annot_prom$Gene[1:10])
for (i in 1:length(topgenes_prom)){
  LIMMAout_subset <- LIMMAout_annot_prom[(LIMMAout_annot_prom$Gene==topgenes_prom[i]) & 
                                           (LIMMAout_annot_prom$adj.P.Val<0.10),]
  if(nrow(LIMMAout_subset)>1){
    print(LIMMAout_subset[sort(LIMMAout_subset$Pos,index.return=T)$ix,])
  }
  
  
}



