#DESeq2 (version 1.12.3)

setwd("/work/agrobacterium/jhonatas/quant_gene_expression/htseq/")

#preparando arquivo do HTSeq
directory = "."
sampleFiles = grep("count.txt",list.files(directory),value = TRUE)
sampleCondition = rep(c("F1","F2","F3","Y1","Y2","Y3"),each=2)
sampleTable = data.frame(sampleName = c("F1","F1_1","F2","F2_1","F3","F3_1","Y1","Y1_1","Y2","Y2_1","Y3","Y3_1"), fileName = sampleFiles, condition = sampleCondition)

#Rodando DESeq2

library(DESeq2)

dds = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~condition)
keepdds = rowSums(counts(dds)) > 1
table(keepdds)
dds = dds[rowSums(counts(dds)) > 1, ] #filtra genes pouco expressos!
dds$condition = droplevels(dds$condition)
dds = DESeq(dds)

#PCA Plot

rld <- rlog(dds, blind=FALSE)
plotPCA(rld)

##Comparações de grupos

### wild_growth_aerobic vs mutant_growth_aerobic
Y1XF1 = results(dds, lfcThreshold = 1, contrast = c("condition","Y1","F1"));
Y1XF1 = Y1XF1[order(Y1XF1$log2FoldChange),] # ordenar por LFC
summary(Y1XF1)

### wild_production_aerobic vs mutant_production_aerobic
Y2XF2 = results(dds, lfcThreshold = 1, contrast = c("condition","Y2","F2"))
Y2XF2 = Y2XF2[order(Y2XF2$log2FoldChange),] # ordenar por LFC
summary(Y2XF2)

### wild_production_anoxic vs mutant_production_anoxic
Y3XF3 = results(dds, lfcThreshold = 1, contrast = c("condition","Y3","F3"))
Y3XF3 = Y3XF3[order(Y3XF3$log2FoldChange),] # ordenar por LFC
summary(Y3XF3)

### wild_growth_aerobic vs wild_production_aerobic
Y1XY2 = results(dds, lfcThreshold = 1, contrast = c("condition","Y1","Y2"))
Y1XY2 = Y1XY2[order(Y1XY2$log2FoldChange),] # ordenar por LFC
summary(Y1XY2)

### wild_production_aerobic vs wild_production_anoxic
Y2XY3 = results(dds, lfcThreshold = 1, contrast = c("condition","Y2","Y3"))
Y2XY3 = Y2XY3[order(Y2XY3$log2FoldChange),] # ordenar por LFC
summary(Y2XY3)

### mutant_growth_aerobic vs mutant_production_aerobic
F1XF2 = results(dds, lfcThreshold = 1, contrast = c("condition","F1","F2"))
F1XF2 = F1XF2[order(F1XF2$log2FoldChange),] # ordenar por LFC
summary(F1XF2)

### mutant_production_aerobic vs mutant_production_anoxic
F2XF3 = results(dds, lfcThreshold = 1, contrast = c("condition","F2","F3"))
F2XF3 = F2XF3[order(F2XF3$log2FoldChange),] # ordenar por LFC
summary(F2XF3)

##Preparando tabela de LFC para Heatmap, ordenar results pelo nome do gene!!!

Y1XY1 = rep_len(0,length.out = length(rownames(counts(dds))))

Y1XY2 = results(dds, lfcThreshold = 1, contrast = c("condition","Y1","Y2"))
Y1XY2 = Y1XY2[order(rownames(Y1XY2)),] #ordenar pelo nome do gene

Y1XY3 = results(dds, lfcThreshold = 1, contrast = c("condition","Y1","Y3"))
Y1XY3 = Y1XY3[order(rownames(Y1XY2)),] #ordenar pelo nome do gene

Y1XF1 = results(dds, lfcThreshold = 1, contrast = c("condition","Y1","F1"))
Y1XF1 = Y1XF1[order(rownames(Y1XY2)),] #ordenar pelo nome do gene

Y1XF2 = results(dds, lfcThreshold = 1, contrast = c("condition","Y1","F2"))
Y1XF2 = Y1XF2[order(rownames(Y1XY2)),] #ordenar pelo nome do gene

Y1XF3 = results(dds, lfcThreshold = 1, contrast = c("condition","Y1","F3"))
Y1XF3 = Y1XF3[order(rownames(Y1XY2)),] #ordenar pelo nome do gene

DE_genes <- read.table("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/DE_genes.txt", quote="\"", comment.char="")
DE_genes = as.character(DE_genes$V1)

lfc_DEG = data.frame(Y1XY1,Y1XY2$log2FoldChange,Y1XY3$log2FoldChange,Y1XF1$log2FoldChange,Y1XF2$log2FoldChange,Y1XF3$log2FoldChange)
colnames(lfc_DEG) = c("Y1vsY1","Y1vsY2","Y1vsY3","Y1vsF1","Y1vsF2","Y1vsF3")
rownames(lfc_DEG) = rownames(Y1XY2)
lfc_DEG = lfc_DEG[DE_genes,]
setwd("./DEG/deseq2/lfc_threshold/")
save(lfc_DEG, file = "lfc_DEG.RData")

#Exportando resultados

setwd("./DEG/deseq2/lfc_threshold/")

Y1XF1_out = subset(Y1XF1, Y1XF1$padj <= 0.05)
write.table(as.data.frame(Y1XF1_out), file = "wild_growth_aerobicXmutant_growth_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

Y2XF2_out = subset(Y2XF2, Y2XF2$padj <= 0.05)
write.table(as.data.frame(Y2XF2_out), file = "wild_production_aerobicXmutant_production_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

Y3XF3_out = subset(Y3XF3, Y3XF3$padj <= 0.05)
write.table(as.data.frame(Y3XF3_out), file = "wild_production_anoxicXmutant_production_anoxic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

Y1XY2_out = subset(Y1XY2, Y1XY2$padj <= 0.05)
write.table(as.data.frame(Y1XY2_out), file = "wild_growth_aerobicXwild_production_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

Y2XY3_out = subset(Y2XY3, Y2XY3$padj <= 0.05)
write.table(as.data.frame(Y2XY3_out), file = "wild_production_aerobicXwild_production_anoxic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

F1XF2_out = subset(F1XF2, F1XF2$padj <= 0.05)
write.table(as.data.frame(F1XF2_out), file = "mutant_growth_aerobicXmutant_production_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

F2XF3_out = subset(F2XF3, F2XF3$padj <= 0.05)
write.table(as.data.frame(F2XF3_out), file = "mutant_production_aerobicXmutant_production_anoxic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")
