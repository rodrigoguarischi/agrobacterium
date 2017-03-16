#clusterProfiler (version 3.0.5)

options(stringsAsFactors = FALSE)

#Instalando clusterProfiler

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
require(DOSE)
require(clusterProfiler)

#Carregando anotações com COG

protein.id_cog <- read.delim( 
  file = "/mnt/work1/agrobacterium/jhonatas/annotation/cog/cdd2cog/results.bkp/protein-id_cog.txt", 
  header=FALSE, 
  sep = "\t",
  col.names = 1:max( count.fields("/mnt/work1/agrobacterium/jhonatas/annotation/cog/cdd2cog/results.bkp/protein-id_cog.txt", sep = "\t") )
  );

cog_stats <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/cog/cdd2cog/results.bkp/cog_stats.txt", header=FALSE)

cog2gene = protein.id_cog[,c(2,1)]
cog2name = cog_stats[,c("V1","V2")]

#Carregando anotações com KEGG

kegg2gene <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/kegg/results/kegg2gene.txt")
kegg2name <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/kegg/results/kegg2name.txt")

#Carregando anotações com GO

go2gene <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/interproscan/results/go2gene.tsv")
go2name <- read.delim("/mnt/work1/agrobacterium/jhonatas/annotation/interproscan/results/go2name.tsv", header=FALSE)

#Carregando lista de genes diferencialmente expressos para rodar análise de enriquecimento com COG 

for( file in list.files(path = "/work/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold", pattern = "*txt", full.names = T) ){
  deg = rownames(read.delim( file = file, header = T))

  #Executando o pacote clusterProfiler para cada "file"
  x = enricher(deg, TERM2GENE = cog2gene, TERM2NAME = cog2name)
  
  #Apresentando os resultados
  
  print( paste0("Arquivo = " , file) );

  #print( table( deg %in% as.character( cog2gene[,2] ) ) )
    
  print( summary(x) )
  
  print(barplot(x))

}

#Carregando lista de genes diferencialmente expressos para rodar análise de enriquecimento com KEGG 

for( file in list.files(path = "/work/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold", pattern = "*txt", full.names = T) ){
  deg = rownames(read.delim( file = file, header = T))
  
  #Executando o pacote clusterProfiler para cada "file"
  x = enricher(deg, TERM2GENE = kegg2gene, TERM2NAME = kegg2name)
  
  #Apresentando os resultados
  
  print( paste0("Arquivo = " , file) );
  
  #print( table( deg %in% as.character( cog2gene[,2] ) ) )
  
  print( summary(x) )
  
  print(barplot(x))
  
}

#Carregando lista de genes diferencialmente expressos para rodar análise de enriquecimento com GO 

for( file in list.files(path = "/work/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold", pattern = "*txt", full.names = T) ){
  deg = rownames(read.delim( file = file, header = T))
  
  #Executando o pacote clusterProfiler para cada "file"
  x = enricher(deg, TERM2GENE = go2gene, TERM2NAME = go2name)
  
  #Apresentando os resultados
  
  print( paste0("Arquivo = " , file) );
  
  #print( table( deg %in% as.character( cog2gene[,2] ) ) )
  
  print( summary(x) )
  
  print(barplot(x))
  
}