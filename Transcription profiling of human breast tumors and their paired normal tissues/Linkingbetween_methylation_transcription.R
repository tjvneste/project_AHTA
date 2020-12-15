#set working directory
setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")

load("Significant_output_annotation_transcription.Rda")# significant_pvalues2
load('methylation_all_pvalues.Rda') # adj_pvalues
load("adjusted_pvalues_annotation_transcription.Rda") # adjusted pvalues
load("Methylation_significant.Rda") # significant_p_values


down_RNA<- read.table('down-2.txt', header=FALSE)
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
# this means higher expression in the tumor samples
neg_logfold<- significant_pvalues2[significant_pvalues2$logFC <= (-1),]
# this means higher expression in the control samples
write.table(pos_logfold$hgnc_symbol,file='transcription_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$hgnc_symbol,file='transcription_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

#writing tables => for all transcription significant genes
write.table(adjusted_pvalues$hgnc_symbol,file='transcription_genes_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
pos_expression <- adjusted_pvalues[adjusted_pvalues$logFC >= 0,] 
length(pos_expression$hgnc_symbol) #3211
write.table(pos_expression$hgnc_symbol,file='transcription_genes_pos_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
neg_expression <- adjusted_pvalues[adjusted_pvalues$logFC <= 0,] 
length(neg_expression$hgnc_symbol) # 1817
write.table(neg_expression$hgnc_symbol,file='transcription_genes_neg_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

#writing tables => for methylation 622 significant genes
write.table(significant_p_values$Gene,file='methylation_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
pos_logfold<- significant_p_values[significant_p_values$logFC >= 2,]
dim(pos_logfold) # 481
# this means higher expression in the tumor samples
neg_logfold<- significant_p_values[significant_p_values$logFC <= (-2),]
dim(neg_logfold)
# this means higher expression in the control samples
write.table(pos_logfold$Gene,file='methylation_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$Gene,file='methylation_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

#writing tables => for methylation 13179 significant genes
write.table(adj_pvalues$Gene,file='methylation_genes_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
pos_logfold<- adj_pvalues[adj_pvalues$logFC > 0,]
dim(pos_logfold) # 8054
# this means higher expression in the tumor samples
neg_logfold<- adj_pvalues[adj_pvalues$logFC < 0,]
dim(neg_logfold) # 5125
# this means higher expression in the control samples
write.table(pos_logfold$Gene,file='methylation_genes_pos_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$Gene,file='methylation_genes_neg_logfold_all.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

sessionInfo()
