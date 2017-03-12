#DESeq2 (version 1.12.3)
#tximport (version 1.0.3)

setwd("/work/agrobacterium/jhonatas/quant_gene_expression/salmon/")

#preparando arquivo do Salmon com Tximport:

library(tximport)
wild_growth_aerobic_0 = file.path("./aln1/salmon_out", "quant.sf")
wild_growth_aerobic_1 = file.path("./aln2/salmon_out", "quant.sf")
wild_production_aerobic_0 = file.path("./aln3/salmon_out", "quant.sf")
wild_production_aerobic_1 = file.path("./aln4/salmon_out", "quant.sf")
wild_production_anoxic_0 = file.path("./aln5/salmon_out", "quant.sf")
wild_production_anoxic_1 = file.path("./aln6/salmon_out", "quant.sf")
mutant_growth_aerobic_0 = file.path("./aln7/salmon_out", "quant.sf")
mutant_growth_aerobic_1 = file.path("./aln8/salmon_out", "quant.sf")
mutant_production_aerobic_0 = file.path("./aln9/salmon_out", "quant.sf")
mutant_production_aerobic_1 = file.path("./aln10/salmon_out", "quant.sf")
mutant_production_anoxic_0 = file.path("./aln11/salmon_out", "quant.sf")
mutant_production_anoxic_1 = file.path("./aln12/salmon_out", "quant.sf")

files = c(wild_growth_aerobic_0,wild_growth_aerobic_1,wild_production_aerobic_0,wild_production_aerobic_1,
          wild_production_anoxic_0,wild_production_anoxic_1,mutant_growth_aerobic_0,mutant_growth_aerobic_1,mutant_production_aerobic_0,
          mutant_production_aerobic_1,mutant_production_anoxic_0,mutant_production_anoxic_1)
names(files) = c("wga","wga_1","wpa","wpa_1","wpn","wpn_1","mga","mga_1","mpa","mpa_1","mpn","mpn_1")

txi = tximport(files, type = "salmon", txOut = TRUE)

##grupos separados
#sampleTable = data.frame(type = factor(rep(c("wild","mutant"),each=6)), condition = factor(c("growth","growth","production","production","production","production","growth","growth","production","production","production","production")), treatment = factor(c("aerobic","aerobic","aerobic","aerobic","anoxic","anoxic","aerobic","aerobic","aerobic","aerobic","anoxic","anoxic")))
#rownames(sampleTable) = colnames(txi$counts)

#dds = DESeqDataSetFromTximport(txi, sampleTable, ~ type + condition + treatment)
#keepdds = rowSums(counts(dds)) > 1
#table(keepdds)
#dds = dds[rowSums(counts(dds)) > 1, ] #filtra genes pouco expressos! não precisou...
#dds = DESeq(dds)

##grupos únicos, estou usando esse!

sampleTable = data.frame(group = factor(rep(c("wga","wpa","wpn","mga","mpa","mpn"),each=2)))
rownames(sampleTable) = colnames(txi$counts)

#Rodando DESeq2

library(DESeq2)

dds = DESeqDataSetFromTximport(txi, sampleTable, ~ group)
keepdds = rowSums(counts(dds)) > 1
table(keepdds)
#dds = dds[rowSums(counts(dds)) > 1, ] #filtra genes pouco expressos! não precisou...
dds = DESeq(dds)

#PCA Plot

rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup = c("condition"), ntop = 5161)

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


#Exportando resultados

setwd("./DEG/")

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
