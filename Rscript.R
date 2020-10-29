setwd("C:/Users/Boris/Documents/School/2de master/Applied high throughput analysis/Practical sessions/PC4")

rnadatafull <- read.table("SalivaryGlandCancer.txt", sep="\t", header=TRUE, quote="")

rnadatafull <- rnadatafull[duplicated(rnadatafull[,5:8])==F,]
rownames(rnadatafull) <- 1:nrow(rnadatafull)
rnadata <- rnadatafull[,9:14]

library("edgeR")
y <-DGEList(rnadata, genes=rnadatafull[1])
y <- calcNormFactors(y)
y$samples$lib.size <- colSums(y$counts)

#cutoff <- 3/(mean(y$samples$lib.size)/100000)
#Keep rows where cpm(y) > cutoff
#y <- keep[y, , ]

pdf("MDS.pdf")
plotMDS(y)
dev.off()

patient <- factor(c('A', 'B', 'C','A', 'B', 'C'))
case <- factor(c('_control', '_control', '_control', '_leasion', '_leasion', '_leasion'))
data.frame(sample=colnames(y), patient, case)
design <- model.matrix(~patient+case)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design, robust=TRUE)
pdf("BCVplot.pdf")
plotBCV(y)
dev.off()

fit <- glmQLFit(y, design)
lrt <- glmQLFTest(fit, coef=2)
topTags(lrt)

summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

BiocManager::install("DESeq2")
library("DESeq2")

cond <- factor(rep(c("control","tumor"),each=3))
ind <- factor(rep(c(1,2,3),2))
colData <- DataFrame(ind,cond)

ddsP <- DESeqDataSetFromMatrix(rnadata,colData,~ind+cond)
dss <- DESeq(ddsP)

png("qc-dispersions.png", 1500, 1500, pointsize=20) #maak een afbeelding, 20 is de tekstgrootte
plotDispEsts(dds, main="dispersion plot")
dev.off()

normcnt <- counts(dds,normalized=T)
rld <- rlogTransformation(dds)


