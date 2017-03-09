#Abrir arquivo com GO IDs

gene2GO <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/interproscan/gene2GO.tsv", header=FALSE)
colnames(gene2GO) = c("gene","go_id")

#Mapeando GO Ids com respectivos termos e ontologia

library(clusterProfiler)

GO2term = go2term(as.character(gene2GO$V2))
GO2Ont = go2ont(as.character(gene2GO$V2))
gene_go_term = merge(gene2GO,GO2term, by = "go_id", all = TRUE)
go_annotation = merge(gene_go_term, GO2Ont, by = "go_id", all = TRUE )

#Removendo arquivos intermediários

rm(gene_go_term)
rm(GO2Ont)
rm(GO2term)
rm(gene2GO)

#Exportando os resultados

write.table(go_annotation, file = "go_annotation.tsv", quote=FALSE, row.names=FALSE, col.names = TRUE, sep = "\t")

##Remover elementos repetidos no shell depois com o comando sort "go_annotation.tsv | uniq > final_go_annotation.tsv" e "mv final_go_annotation.tsv go_annotation.tsv"