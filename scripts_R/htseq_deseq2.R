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
wgaXmga = results(dds, lfcThreshold = 1, contrast = c("condition","wga","mga"));
wgaXmga = wgaXmga[order(wgaXmga$log2FoldChange),]
summary(wgaXmga)

### wild_production_aerobic vs mutant_production_aerobic
wpaXmpa = results(dds, lfcThreshold = 1, contrast = c("condition","wpa","mpa"))
wpaXmpa = wpaXmpa[order(wpaXmpa$log2FoldChange),]
summary(wpaXmpa)

### wild_production_anoxic vs mutant_production_anoxic
wpnXmpn = results(dds, lfcThreshold = 1, contrast = c("condition","wpn","mpn"))
wpnXmpn = wpnXmpn[order(wpnXmpn$log2FoldChange),]
summary(wpnXmpn)

### wild_growth_aerobic vs wild_production_aerobic
wgaXwpa = results(dds, lfcThreshold = 1, contrast = c("condition","wga","wpa"))
wgaXwpa = wgaXwpa[order(wgaXwpa$log2FoldChange),]
summary(wgaXwpa)

### wild_production_aerobic vs wild_production_anoxic
wpaXwpn = results(dds, lfcThreshold = 1, contrast = c("condition","wpa","wpn"))
wpaXwpn = wpaXwpn[order(wpaXwpn$log2FoldChange),]
summary(wpaXwpn)

### mutant_growth_aerobic vs mutant_production_aerobic
mgaXmpa = results(dds, lfcThreshold = 1, contrast = c("condition","mga","mpa"))
mgaXmpa = mgaXmpa[order(mgaXmpa$log2FoldChange),]
summary(mgaXmpa)

### mutant_production_aerobic vs mutant_production_anoxic
mpaXmpn = results(dds, lfcThreshold = 1, contrast = c("condition","mpa","mpn"))
mpaXmpp = mpaXmpn[order(mpaXmpn$log2FoldChange),]
summary(mpaXmpn)

#Heatmap (version 1.0.8)
library(pheatmap)
tab_contagem = as.data.frame(counts(dds)) #tabela com as contagens
DEgenes = read.table("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/DEgenes.txt", quote="\"", comment.char="")
DEgenes = as.character(DEgenes$V1) #lista de todos os genes DE em todas as condições
DE_counts = tab_contagem[DEgenes,] #tabela com as contagens para os genes DE
log2_DEcounts = log2(DE_counts+1) #transformando contagens pelo Log2
sampleDist = dist(t(log2_DEcounts)) #Calculando distância euclidiana entre as amostras
sampleDistMatrix = as.matrix(sampleDist) #Gerando matriz de distâncias

annotation = data.frame(Condition = rep(c("F1","F2","F3","Y1","Y2","Y3"),each=2))
rownames(annotation) = c("F1","F1_1","F2","F2_1","F3","F3_1","Y1","Y1_1","Y2","Y2_1","Y3","Y3_1")
ann_color = list(Condition = c(F1 ="blue", F2 = "yellow", F3 = "magenta", Y1 = "brown", Y2 = "darkcyan", Y3 = "lavender"))

library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(9,"Reds")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDist, clustering_distance_cols = sampleDist,
         color = color, legend = TRUE, annotation_row = annotation, annotation_col = annotation, annotation_colors = ann_color)

##install.packages("BoutrosLab.plotting.general", repos = "http://bpg.oicr.on.ca") ##Instalando BoutrosLab package

#Exportando resultados

setwd("./DEG/deseq2/lfc_threshold/")

wgaXmga_out = subset(wgaXmga, wgaXmga$padj <= 0.05)
write.table(as.data.frame(wgaXmga_out), file = "wild_growth_aerobicXmutant_growth_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

wpaXmpa_out = subset(wpaXmpa, wpaXmpa$padj <= 0.05)
write.table(as.data.frame(wpaXmpa_out), file = "wild_production_aerobicXmutant_production_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

wpnXmpn_out = subset(wpnXmpn, wpnXmpn$padj <= 0.05)
write.table(as.data.frame(wpnXmpn_out), file = "wild_production_anoxicXmutant_production_anoxic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

wgaXwpa_out = subset(wgaXwpa, wgaXwpa$padj <= 0.05)
write.table(as.data.frame(wgaXwpa_out), file = "wild_growth_aerobicXwild_production_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

wpaXwpn_out = subset(wpaXwpn, wpaXwpn$padj <= 0.05)
write.table(as.data.frame(wpaXwpn_out), file = "wild_production_aerobicXwild_production_anoxic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

mgaXmpa_out = subset(mgaXmpa, mgaXmpa$padj <= 0.05)
write.table(as.data.frame(mgaXmpa_out), file = "mutant_growth_aerobicXmutant_production_aerobic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")

mpaXmpn_out = subset(mpaXmpn, mpaXmpn$padj <= 0.05)
write.table(as.data.frame(mpaXmpn_out), file = "mutant_production_aerobicXmutant_production_anoxic.txt", sep = "\t", row.names = TRUE, col.names = TRUE, na = "NA", dec = ".")
