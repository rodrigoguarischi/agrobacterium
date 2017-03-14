#EdgeR (version 3.16.5)
setwd("/work/agrobacterium/jhonatas/quant_gene_expression/htseq/")

library(edgeR)

#preparando arquivo do HTSeq
directory = "."
files = grep("count.txt",list.files(directory),value = TRUE)
grupos = rep(c("mga","mpa","mpn","wga","wpa","wpn"),each=2)
x = readDGE(files, columns = c(1,2), group = grupos)

#Filtrando genes pouco expressos
keep = rowSums(x$counts) > 1
summary(keep)
x = x[keep, , keep.lib.sizes=FALSE]
x$samples$lib.size = colSums(x$counts)

#Renomeando tabelas
rownames(x$samples) = c("mga_0","mga_1","mpa_0","mpa_1","mpn_0","mpn_1","wga_0","wga_1","wpa_0","wpa_1","wpn_0","wpn_1")
colnames(x$counts) = c("mga_0","mga_1","mpa_0","mpa_1","mpn_0","mpn_1","wga_0","wga_1","wpa_0","wpa_1","wpn_0","wpn_1")

#Normaliza por TMM
x = calcNormFactors(x)

#Plot PCA
plotMDS(x)

#Definindo dispersão para réplicas técnicas conforme desenvolvedores do EdgeR

bcv = 0.01

#Comparação entre os grupos

##wild_growth_aerobic vs mutant_growth_aerobic
wgaXmga = exactTest(x, pair = c("wga","mga"), dispersion = bcv^2)
wgaXmga_names = as.data.frame(decideTestsDGE(wgaXmga, lfc = 1))
rownames(wgaXmga_names) = rownames(wgaXmga$table)
wgaXmga_filter = (wgaXmga_names$V1 > 0 | wgaXmga_names$V1 < 0)
wgaXmga_tags = as.data.frame(topTags(wgaXmga, n=Inf))
wgaXmga_tags = wgaXmga_tags[order(rownames(wgaXmga_tags)),]
wgaXmga_out = wgaXmga_tags[wgaXmga_filter,]
wgaXmga_out = wgaXmga_out[order(wgaXmga_out$logFC),]
summary(decideTestsDGE(wgaXmga, lfc = 1))

##wild_production_aerobic vs mutant_production_aerobic
wpaXmpa = exactTest(x, pair = c("wpa","mpa"), dispersion = bcv^2)
wpaXmpa_names = as.data.frame(decideTestsDGE(wpaXmpa, lfc = 1))
rownames(wpaXmpa_names) = rownames(wpaXmpa$table)
wpaXmpa_filter = (wpaXmpa_names$V1 > 0 | wpaXmpa_names$V1 < 0)
wpaXmpa_tags = as.data.frame(topTags(wpaXmpa, n=Inf))
wpaXmpa_tags = wpaXmpa_tags[order(rownames(wpaXmpa_tags)),]
wpaXmpa_out = wpaXmpa_tags[wpaXmpa_filter,]
wpaXmpa_out = wpaXmpa_out[order(wpaXmpa_out$logFC),]
summary(decideTestsDGE(wpaXmpa, lfc = 1))

##wild_production_anoxic vs mutant_production_anoxic
wpnXmpn = exactTest(x, pair = c("wpn","mpn"), dispersion = bcv^2)
wpnXmpn_names = as.data.frame(decideTestsDGE(wpnXmpn, lfc = 1))
rownames(wpnXmpn_names) = rownames(wpnXmpn$table)
wpnXmpn_filter = (wpnXmpn_names$V1 > 0 | wpnXmpn_names$V1 < 0)
wpnXmpn_tags = as.data.frame(topTags(wpnXmpn, n=Inf))
wpnXmpn_tags = wpnXmpn_tags[order(rownames(wpnXmpn_tags)),]
wpnXmpn_out = wpnXmpn_tags[wpnXmpn_filter,]
wpnXmpn_out = wpnXmpn_out[order(wpnXmpn_out$logFC),]
summary(decideTestsDGE(wpnXmpn, lfc = 1))

## wild_growth_aerobic vs wild_production_aerobic
wgaXwpa = exactTest(x, pair = c("wga","wpa"), dispersion = bcv^2)
wgaXwpa_names = as.data.frame(decideTestsDGE(wgaXwpa, lfc = 1))
rownames(wgaXwpa_names) = rownames(wgaXwpa$table)
wgaXwpa_filter = (wgaXwpa_names$V1 > 0 | wgaXwpa_names$V1 < 0)
wgaXwpa_tags = as.data.frame(topTags(wgaXwpa, n=Inf))
wgaXwpa_tags = wgaXwpa_tags[order(rownames(wgaXwpa_tags)),]
wgaXwpa_out = wgaXwpa_tags[wgaXwpa_filter,]
wgaXwpa_out = wgaXwpa_out[order(wgaXwpa_out$logFC),]
summary(decideTestsDGE(wgaXwpa, lfc = 1))

## wild_production_aerobic vs wild_production_anoxic
wpaXwpn = exactTest(x, pair = c("wpa","wpn"), dispersion = bcv^2)
wpaXwpn_names = as.data.frame(decideTestsDGE(wpaXwpn, lfc = 1))
rownames(wpaXwpn_names) = rownames(wpaXwpn$table)
wpaXwpn_filter = (wpaXwpn_names$V1 > 0 | wpaXwpn_names$V1 < 0)
wpaXwpn_tags = as.data.frame(topTags(wpaXwpn, n=Inf))
wpaXwpn_tags = wpaXwpn_tags[order(rownames(wpaXwpn_tags)),]
wpaXwpn_out = wpaXwpn_tags[wpaXwpn_filter,]
wpaXwpn_out = wpaXwpn_out[order(wpaXwpn_out$logFC),]
summary(decideTestsDGE(wpaXwpn, lfc = 1))

## mutant_growth_aerobic vs mutant_production_aerobic
mgaXmpa = exactTest(x, pair = c("mga","mpa"), dispersion = bcv^2)
mgaXmpa_names = as.data.frame(decideTestsDGE(mgaXmpa, lfc = 1))
rownames(mgaXmpa_names) = rownames(mgaXmpa$table)
mgaXmpa_filter = (mgaXmpa_names$V1 > 0 | mgaXmpa_names$V1 < 0)
mgaXmpa_tags = as.data.frame(topTags(mgaXmpa, n=Inf))
mgaXmpa_tags = mgaXmpa_tags[order(rownames(mgaXmpa_tags)),]
mgaXmpa_out = mgaXmpa_tags[mgaXmpa_filter,]
mgaXmpa_out = mgaXmpa_out[order(mgaXmpa_out$logFC),]
summary(decideTestsDGE(mgaXmpa, lfc = 1))

## mutant_production_aerobic vs mutant_production_anoxic
mpaXmpn = exactTest(x, pair = c("mpa","mpn"), dispersion = bcv^2)
mpaXmpn_names = as.data.frame(decideTestsDGE(mpaXmpn, lfc = 1))
rownames(mpaXmpn_names) = rownames(mpaXmpn$table)
mpaXmpn_filter = (mpaXmpn_names$V1 > 0 | mpaXmpn_names$V1 < 0)
mpaXmpn_tags = as.data.frame(topTags(mpaXmpn, n=Inf))
mpaXmpn_tags = mpaXmpn_tags[order(rownames(mpaXmpn_tags)),]
mpaXmpn_out = mpaXmpn_tags[mpaXmpn_filter,]
mpaXmpn_out = mpaXmpn_out[order(mpaXmpn_out$logFC),]
summary(decideTestsDGE(mpaXmpn, lfc = 1))

#Exportando resultados
setwd("./DEG/edger/")

write.table(as.data.frame(wgaXmga_out), file = "wild_growth_aerobic_X_mutant_growth_aerobic.txt", row.names = TRUE, col.names =TRUE, quote = FALSE, sep = "\t", na = "NA", dec = ".")
write.table(as.data.frame(wpaXmpa_out), file = "wild_production_aerobic_X_mutant_production_aerobic.txt", row.names = TRUE, col.names =TRUE, quote = FALSE, sep = "\t", na = "NA", dec = ".")
write.table(as.data.frame(wpnXmpn_out), file = "wild_production_anoxic_X_mutant_production_anoxic.txt", row.names = TRUE, col.names =TRUE, quote = FALSE, sep = "\t", na = "NA", dec = ".")
write.table(as.data.frame(wgaXwpa_out), file = "wild_growth_aerobic_X_wild_production_aerobic.txt", row.names = TRUE, col.names =TRUE, quote = FALSE, sep = "\t", na = "NA", dec = ".")
write.table(as.data.frame(wpaXwpn_out), file = "wild_production_aerobic_X_wild_production_anoxic.txt", row.names = TRUE, col.names =TRUE, quote = FALSE, sep = "\t", na = "NA", dec = ".")
write.table(as.data.frame(mgaXmpa_out), file = "mutant_growth_aerobic_X_mutant_production_aerobic.txt", row.names = TRUE, col.names =TRUE, quote = FALSE, sep = "\t", na = "NA", dec = ".")
write.table(as.data.frame(mpaXmpn_out), file = "mutant_production_aerobic_X_mutant_production_anoxic.txt", row.names = TRUE, col.names =TRUE, quote = FALSE, sep = "\t", na = "NA", dec = ".")