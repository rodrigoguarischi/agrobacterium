#Instalando KEGGREST

source("https://bioconductor.org/biocLite.R")
biocLite("KEGGREST")

#Chamando pacote KEGGREST

library(KEGGREST)

#Carregando dados
KO_list = read.table("/mnt/work1/agrobacterium/jhonatas/annotation/kegg/KO_list.txt", quote="\"", comment.char="")
kos.ids = as.character(KO_list[,1])

#Rodando KEGGREST

query <- NULL;

for ( i in seq(1, length( kos.ids ), 10 ) ){
  start <- i;
  end <- ifelse(i+9 <= length( kos.ids ), i+9, length( kos.ids ) );
  print( paste0( "querying from " , start , " to ", end ) );
  query <- c(query, keggGet( kos.ids[start:end] ) )
}

#Extraindo definições para os KEGG-IDs e vias metabólicas envolvidas

ENTRY = NULL;
DEFINITION = NULL;
PATHWAY = NULL;

for ( i in 1:length(kos.ids)){
    ENTRY = c(ENTRY, query[[i]]$ENTRY);
    DEFINITION = c(DEFINITION, query[[i]]$DEFINITION);
    PATHWAY = c(PATHWAY, c(kos.ids[i],query[[i]]$PATHWAY))
}

paths_description = stack(PATHWAY)
keggID2keggTerm = data.frame(ENTRY,DEFINITION)

#Formatar "Path_description.tsv" com o script "parse_pathwaysKEGG.pl" e importar a saída para o R

write.table(paths_description, file = "Path_description.tsv", quote=FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")
keggID2pathways <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/kegg/keggID2pathways.tsv", header=FALSE)
colnames(keggID2pathways) = c("ENTRY", "PATHWAYS", "PATHWAYS-ID")

#Juntar os arquivos de saída em uma única tabela

kegg_annotation = merge(keggID2keggTerm,keggID2pathways, by="ENTRY", all = TRUE)

#Mapear genes com a tabela de anotação do KEGG

genes2keggID_only_with_keggID <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/kegg/genes2keggID_only_with_keggID.txt", header=FALSE)
colnames(genes2keggID_only_with_keggID) = c("GENE","ENTRY")
genes_kegg_annotation = merge(genes2keggID_only_with_keggID,kegg_annotation, by="ENTRY", all = TRUE)
colnames(genes_kegg_annotation) = c("KEGG-ID","GENE","KEGG-TERM","PATHWAY","PATHWAY-ID")

#Removendo arquivos intermediários
rm(paths_description)
rm(genes2keggID_only_with_keggID)
rm(KO_list)
rm(kegg_annotation)
rm(keggID2pathways)
rm(keggID2keggTerm)

#Exportando os resultados

write.table(genes_kegg_annotation, file = "genes_kegg_annotation.tsv", quote=FALSE, row.names=FALSE, col.names = TRUE, sep = "\t")
