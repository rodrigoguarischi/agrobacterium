#Instalando clusterProfiler

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
require(DOSE)
require(clusterProfiler)

#Carregando anotações com COG

protein.id_cog <- read.delim("/mnt/work1/agrobacterium/jhonatas/cdd2cog/cdd2cog/results.bkp/protein-id_cog.txt", header=FALSE)
cog_stats <- read.delim("/mnt/work1/agrobacterium/jhonatas/cdd2cog/cdd2cog/results.bkp/cog_stats.txt", header=FALSE)
cog2gene = protein.id_cog[,c("V2","V1")]
cog2name = cog_stats[,c("V1","V2")]

#Carregando lista de genes diferencialmente expressos

deg = as.vector(list_genes$V1)

#Executando o pacote clusterProfiler

x = enricher(deg, TERM2GENE = cog2gene, TERM2NAME = cog2name, qvalueCutoff = 0.05 )

#Apresentando os resultados

summary(x)
barplot(x)