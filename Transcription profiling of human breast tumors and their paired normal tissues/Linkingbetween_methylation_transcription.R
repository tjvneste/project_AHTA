load("Significant_output_annotation_transcription.Rda")# object noemt significant_pvalues
load("Significant_output_annotation_methylation.Rda") #LIMMAout_annot_2

dim(LIMMAout_annot_2)
dim(significant_pvalues)
# filtering op limma
LIMMAout_annot <- LIMMAout_annot_2[LIMMAout_annot_2$adj.P.Val<=0.0615913,]
dim(LIMMAout_annot)
dim(significant_pvalues)

head(LIMMAout_annot)
head(significant_pvalues)

sum(LIMMAout_annot$Gene%in%significant_pvalues$hgnc_symbol) # 321 genen die in allebei voorkomen
