setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")

load("Significant_output_annotation_transcription.Rda")# significant_pvalues2
load("Methylation_significant.Rda") # significant_p_values

dim(significant_pvalues2) # 39 11
dim(significant_p_values) #250 8

head(significant_pvalues2) # transcription 
head(significant_p_values) # methylation

sum(significant_pvalues2$hgnc_symbol%in%significant_p_values$Gene) # 0 common genes 

write.table(significant_pvalues2$hgnc_symbol,file='transcription_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

pos_logfold<- significant_pvalues2[significant_pvalues2$logFC >= 1,]
# this means higher expression in the tumor samples
neg_logfold<- significant_pvalues2[significant_pvalues2$logFC <= (-1),]
# this means higher expression in the control samples

write.table(pos_logfold$hgnc_symbol,file='transcription_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$hgnc_symbol,file='transcription_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
## Interpretation results

sessionInfo()
