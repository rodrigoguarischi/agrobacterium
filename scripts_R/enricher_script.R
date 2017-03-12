#clusterProfiler (version 3.0.5)

options(stringsAsFactors = FALSE)

#Instalando clusterProfiler

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
require(DOSE)
require(clusterProfiler)

#Carregando anotações com COG

protein.id_cog <- read.delim( 
  file = "/mnt/work1/agrobacterium/jhonatas/annotation/cdd2cog/cdd2cog/results.bkp/protein-id_cog.txt", 
  header=FALSE, 
  sep = "\t",
  col.names = 1:max( count.fields("/mnt/work1/agrobacterium/jhonatas/annotation/cdd2cog/cdd2cog/results.bkp/protein-id_cog.txt", sep = "\t") )
  );

cog_stats <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/cdd2cog/cdd2cog/results.bkp/cog_stats.txt", header=FALSE)

cog2gene = protein.id_cog[,c(2,1)]
cog2name = cog_stats[,c("V1","V2")]

#Carregando anotações com KEGG

kegg2gene <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/kegg/results/kegg2gene.txt")
kegg2name <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/kegg/results/kegg2name.txt")

#Carregando lista de genes diferencialmente expressos

for( file in list.files(path = "new", pattern = "new_*", full.names = T) ){
  deg = rownames(read.delim( file = file, header = T))

  #Executando o pacote clusterProfiler para cada "file"
  x = enricher( deg, TERM2GENE = cog2gene, TERM2NAME = cog2name, minGSSize = 2, pAdjustMethod = "none")
  
  #Apresentando os resultados
  
  print( paste0("Arquivo = " , file) );

  print( table( deg %in% as.character( cog2gene[,2] ) ) )
    
  print( summary(x) )
  
  print( " " )
  print( " " )
  print( " " )

}

