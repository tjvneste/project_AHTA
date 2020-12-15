#set working directory
setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")

library(goseq)
library(tidyverse)
library(ggplot2)

BiocManager::install("org.Mm.eg.db")

load("Significant_output_annotation_transcription.Rda")# significant_pvalues2
load('methylation_all_pvalues.Rda') # adj_pvalues
load("adjusted_pvalues_annotation_transcription.Rda") # adjusted pvalues
load("Methylation_significant.Rda") # significant_p_values


down_RNA<- read.table('down-2.txt', header=FALSE,quote= "")
up_RNA<- read.table('up-2.txt', header=FALSE)
 
dim(down_RNA) # 184
dim(up_RNA)  #239

total_RNA <- rbind(down_RNA,up_RNA)
head(total_RNA)
write.table(total_RNA,file='RNA_seq_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
dim(total_RNA)

dim(significant_pvalues2) # 39 11
# significant adusted p values(0.05) for transcription and log fold changes of > 1
dim(adjusted_pvalues) # 5028
# all significant adjusted p values(0.05) for transcription
dim(significant_p_values) #622 8
# signifcant adjusted p values(0.010) for methylation and log fold changes of > 2
dim(adj_pvalues) # 13179
# all significant adjusted p values(0.10) for methylation 

head(up_RNA)
head(adj_pvalues)
head(adjusted_pvalues) 
head(significant_pvalues2) # transcription 
head(significant_p_values) # methylation

# writing tables => for transcription 39 significant genes 
write.table(significant_pvalues2$hgnc_symbol,file='transcription_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
pos_logfold<- significant_pvalues2[significant_pvalues2$logFC >= 1,]
dim(pos_logfold) # 9 11
# this means higher expression in the tumor samples
neg_logfold<- significant_pvalues2[significant_pvalues2$logFC <= (-1),]
dim(neg_logfold) #30 11
# this means higher expression in the control samples
write.table(pos_logfold$hgnc_symbol,file='transcription_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$hgnc_symbol,file='transcription_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

#writing tables => for all transcription significant genes
write.table(adjusted_pvalues$hgnc_symbol,file='transcription_genes_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
length(unique(adjusted_pvalues$hgnc_symbol)) # 4069
pos_expression <- adjusted_pvalues[adjusted_pvalues$logFC >= 0,] 
length(unique(pos_expression$hgnc_symbol)) #2621
write.table(pos_expression$hgnc_symbol,file='transcription_genes_pos_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
neg_expression <- adjusted_pvalues[adjusted_pvalues$logFC <= 0,] 
length(unique(neg_expression$hgnc_symbol)) # 1495
write.table(neg_expression$hgnc_symbol,file='transcription_genes_neg_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

#writing tables => for methylation 622 significant genes
write.table(significant_p_values$Gene,file='methylation_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
length(unique(significant_p_values$Gene)) #348
dim(significant_p_values)
pos_logfold<- significant_p_values[significant_p_values$logFC >= 2,]
length(unique(pos_logfold$Gene)) # 250
# this means higher expression in the tumor samples
neg_logfold<- significant_p_values[significant_p_values$logFC <= (-2),]
length(unique(neg_logfold$Gene)) # 104
# this means higher expression in the control samples
write.table(pos_logfold$Gene,file='methylation_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$Gene,file='methylation_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

#writing tables => for methylation 13179 significant genes
write.table(adj_pvalues$Gene,file='methylation_genes_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
length(unique(adj_pvalues$Gene)) # 5917
pos_logfold<- adj_pvalues[adj_pvalues$logFC > 0,]
length(unique(pos_logfold$Gene)) # 3862
# this means higher expression in the tumor samples
neg_logfold<- adj_pvalues[adj_pvalues$logFC < 0,]
length(unique(neg_logfold$Gene)) # 2744
# this means higher expression in the control samples
write.table(pos_logfold$Gene,file='methylation_genes_pos_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$Gene,file='methylation_genes_neg_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)


# Is methylation in promotor region correlated with down regulation of expression? 
## Select probes in promoter regions for whole dataset
LIMMAout_annot_prom_all <- adj_pvalues[grepl("TSS",adj_pvalues$Feature) | (adj_pvalues$Feature=="1stExon"),] #TSS = transcription start site
length(unique(LIMMAout_annot_prom_all$Gene)) # 2100 
head(LIMMAout_annot_prom_all)


positive_methylation_promotor <- LIMMAout_annot_prom_all[LIMMAout_annot_prom_all$logFC > 0,]
negative_methylation_promotor <- LIMMAout_annot_prom_all[LIMMAout_annot_prom_all$logFC < 0,]

write.table(positive_methylation_promotor$Gene,file='methylation_promotor_pos.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(negative_methylation_promotor$Gene,file='methylation_promotor_neg.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
length(unique(positive_methylation_promotor$Gene)) # 1124 
length(unique(negative_methylation_promotor$Gene)) # 999 

# => more methylation in tumor promotor
sum(unique(positive_methylation_promotor$Gene)%in%unique(neg_expression$hgnc_symbol)) # 125
sum(unique(positive_methylation_promotor$Gene)%in%unique(pos_expression$hgnc_symbol)) # 63

# => more methylation in normal promotor
sum(unique(negative_methylation_promotor$Gene)%in%unique(neg_expression$hgnc_symbol)) # 49
sum(unique(negative_methylation_promotor$Gene)%in%unique(pos_expression$hgnc_symbol)) # 119


# Gene Set Enrichment Analysis with the overlapping genes between microarray and the methylation
length(unique(adj_pvalues$Gene)) # 5917 methylation

length(unique(pos_expression$hgnc_symbol)) #2621
length(unique(neg_expression$hgnc_symbol)) # 1495

sum(unique(pos_expression$hgnc_symbol)%in%unique(adj_pvalues$Gene)) #690 
sum(unique(neg_expression$hgnc_symbol)%in%unique(adj_pvalues$Gene)) # 533

list_pos_genes <- unique(pos_expression$hgnc_symbol)[unique(pos_expression$hgnc_symbol)%in%unique(adj_pvalues$Gene)]
list_neg_genes <- unique(neg_expression$hgnc_symbol)[unique(neg_expression$hgnc_symbol)%in%unique(adj_pvalues$Gene)]

write.table(list_pos_genes,file='list_pos_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE) # more expression in tumor
write.table(list_neg_genes,file='list_neg_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE) # less expression in tumor

sessionInfo()
