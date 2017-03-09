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

#Encontrando a hierarquia das vias metabólicas

brites = NULL

for ( i in seq(1, length( kos.ids ), 10 ) ){
  start <- i;
  end <- ifelse(i+9 <= length( kos.ids ), i+9, length( kos.ids ) );
  print( paste0( "querying from " , start , " to ", end ) );
  brites <- c(brites, keggLink("brite", kos.ids[start:end] ) )
}

br = keggList("brite")
tab_br = stack(br)
tab_br = tab_br[c(-1,-2,-3,-4,-5,-6,-102,-103,-104,-105),] #remove BRITES não informativos
tab_brites = stack(brites)

colnames(tab_br) = c("BRITE","ID-BRITE")
colnames(tab_brites) = c("ID-BRITE", "ID-KEGG")
Tab_Path_Brites = merge(tab_brites, tab_br, by="ID-BRITE", all = TRUE)
Tab_Path_Brites = na.omit(Tab_Path_Brites)
rm(tab_br)
rm(tab_brites)

#Exportando os resultados

write.table(keggID2keggTerm, file = "keggID2keggTerm.tsv", quote=FALSE, row.names=FALSE, col.names = TRUE, sep = "\t")
write.table(Tab_Path_Brites, file = "Tab_Path_Brites.tsv", quote=FALSE, row.names=FALSE, col.names = TRUE, sep = "\t")
write.table(paths_description, file = "Path_description.tsv", quote=FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")
