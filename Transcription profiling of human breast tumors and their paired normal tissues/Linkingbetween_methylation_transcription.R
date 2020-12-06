setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")

load("Significant_output_annotation_transcription.Rda")# significant_pvalues2
load('methylation_all_pvalues.Rda') # adj_pvalues
load("adjusted_pvalues_annotation_transcription.Rda") # adjusted pvalues
load("Methylation_significant.Rda") # significant_p_values

dim(significant_pvalues2) # 39 11
dim(significant_p_values) #250 8
dim(adjusted_pvalues) # 5028
dim(adj_pvalues) # 1304


head(adj_pvalues)
head(adjusted_pvalues) 
head(significant_pvalues2) # transcription 
head(significant_p_values) # methylation

sum(significant_pvalues2$hgnc_symbol%in%significant_p_values$Gene) # 0 common genes 
sum(significant_p_values$Gene%in%adjusted_pvalues$hgnc_symbol) # 56 common genes 
sum(adj_pvalues$Gene%in%adjusted_pvalues$hgnc_symbol) # 323 common genes

write.table(significant_pvalues2$hgnc_symbol,file='transcription_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

pos_logfold<- significant_pvalues2[significant_pvalues2$logFC >= 1,]
# this means higher expression in the tumor samples
neg_logfold<- significant_pvalues2[significant_pvalues2$logFC <= (-1),]
# this means higher expression in the control samples

write.table(pos_logfold$hgnc_symbol,file='transcription_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$hgnc_symbol,file='transcription_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
## Interpretation results

sessionInfo()
