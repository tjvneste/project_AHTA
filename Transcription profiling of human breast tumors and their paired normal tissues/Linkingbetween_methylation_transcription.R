setwd("~/Documents/Bioinformatics/Applied high-throughput analysis/project_AHTA/Transcription profiling of human breast tumors and their paired normal tissues")

load("Significant_output_annotation_transcription.Rda")# significant_pvalues2
load('methylation_all_pvalues.Rda') # adj_pvalues
load("adjusted_pvalues_annotation_transcription.Rda") # adjusted pvalues
load("Methylation_significant.Rda") # significant_p_values


down_RNA<- read.table('down.txt', header=FALSE)
up_RNA<- read.table('up.txt', header=FALSE)
 
dim(down_RNA) # 497
dim(up_RNA)  #536

total_RNA <- rbind(down_RNA,up_RNA)
dim(total_RNA)

dim(significant_pvalues2) # 39 11
# dit zijn de significante p values voor transcription and log fold changes van > 1
dim(significant_p_values) #250 8
# dit zijn de signifcante p values voor de methylation and log fold changes van > 2
dim(adjusted_pvalues) # 5028
# alle significante pvalues voor transcription
dim(adj_pvalues) # 1304
# alle significante pvalues voor methylation 

head(up_RNA)
head(adj_pvalues)
head(adjusted_pvalues) 
head(significant_pvalues2) # transcription 
head(significant_p_values) # methylation

write.table(significant_pvalues2$hgnc_symbol,file='transcription_genes.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

pos_logfold<- significant_pvalues2[significant_pvalues2$logFC >= 1,]
# this means higher expression in the tumor samples
neg_logfold<- significant_pvalues2[significant_pvalues2$logFC <= (-1),]
# this means higher expression in the control samples

write.table(pos_logfold$hgnc_symbol,file='transcription_genes_pos_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(neg_logfold$hgnc_symbol,file='transcription_genes_neg_logfold.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

#expression <=> methylation
sum(significant_pvalues2$hgnc_symbol%in%significant_p_values$Gene) # 0 common genes met log foldchange filter
sum(significant_p_values$Gene%in%adjusted_pvalues$hgnc_symbol) # 56 common genes enkel filter op methylation
sum(adj_pvalues$Gene%in%adjusted_pvalues$hgnc_symbol) # 323 common genes geen filter op beide 
#expression <=> RNA
sum(total_RNA$V1 %in% significant_pvalues2$hgnc_symbol) # slechts 31 genes
sum(total_RNA$V1 %in% adjusted_pvalues$hgnc_symbol) # 238 genes
#methylation <=> RNA
sum(total_RNA$V1 %in% significant_p_values$Gene) # slechts 5 genes 
sum(total_RNA$V1 %in% adj_pvalues$Gene) # slechts 31 genes 

# no filter
# postive log fold change
dim(unique(up_RNA))
write.table(unique(up_RNA),file='postive_logfoldchange_RNA.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
#expression
pos_expression <- adjusted_pvalues[adjusted_pvalues$logFC >= 0,] # 3211
length(unique(pos_expression$hgnc_symbol))
write.table(unique(pos_expression$hgnc_symbol),file='postive_logfoldchange_expression.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
#methylation
pos_methylation <- adj_pvalues[adj_pvalues$logFC >= 0,] # 882
pos_methylation_promoter <- pos_methylation[grepl("TSS",pos_methylation$Feature) | (pos_methylation$Feature=="1stExon"),] #TSS = transcription start site
length(unique(pos_methylation_promoter$Gene)) # 148
write.table(unique(pos_methylation$Gene),file='postive_logfoldchange_methylation.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(unique(pos_methylation_promoter$Gene),file='postive_logfoldchange_methylation_promoter.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

# negative log fold change
down_RNA
write.table(unique(down_RNA),file='negative_logfoldchange_RNA.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
#expression
neg_expression <- adjusted_pvalues[adjusted_pvalues$logFC <= 0,] # 1817
length(unique(neg_expression$hgnc_symbol)) # 1495
write.table(unique(neg_expression$hgnc_symbol),file='negative_logfoldchange_expression.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
#methylation
neg_methylation <- adj_pvalues[adj_pvalues$logFC <= 0,] # 422
length(unique(neg_methylation$Gene)) # 296
neg_methylation_promoter <- neg_methylation[grepl("TSS",neg_methylation$Feature) | (neg_methylation$Feature=="1stExon"),] #TSS = transcription start site
length(unique(neg_methylation_promoter$Gene)) # 92
write.table(unique(neg_methylation$Gene),file='negative_logfoldchange_methylation.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(unique(neg_methylation_promoter$Gene),file='negative_logfoldchange_methylation_promoter.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)


sessionInfo()
