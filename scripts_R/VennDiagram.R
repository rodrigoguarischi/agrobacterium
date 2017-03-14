#VennDiagram (version 1.6.17)
library(VennDiagram)

#Carregando arquivos do EdgeR
mutant_growth_aerobic_X_mutant_production_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/mutant_growth_aerobic_X_mutant_production_aerobic.txt")
mutant_production_aerobic_X_mutant_production_anoxic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/mutant_production_aerobic_X_mutant_production_anoxic.txt")
wild_growth_aerobic_X_mutant_growth_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/wild_growth_aerobic_X_mutant_growth_aerobic.txt")
wild_growth_aerobic_X_wild_production_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/wild_growth_aerobic_X_wild_production_aerobic.txt")
wild_production_aerobic_X_mutant_production_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/wild_production_aerobic_X_mutant_production_aerobic.txt")
wild_production_aerobic_X_wild_production_anoxic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/wild_production_aerobic_X_wild_production_anoxic.txt")
wild_production_anoxic_X_mutant_production_anoxic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/wild_production_anoxic_X_mutant_production_anoxic.txt")

#Gerando os diagramas de Venn

setwd("/work/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/edger/vennDiagram/")

##wild_growth_aerobic vs mutant_growth_aerobic, wild_production_aerobic vs mutant_production_aerobic, 
    ## wild_production_anoxic vs mutant_production_anoxic.

venn.diagram(x=list("wgaXmga" = rownames(wild_growth_aerobic_X_mutant_growth_aerobic),
                    "wpaXmpa" = rownames(wild_production_aerobic_X_mutant_production_aerobic),
                    "wpnXmpn" = rownames(wild_production_anoxic_X_mutant_production_anoxic)), filename = "wild_X_mutant.png")

##wild_growth_aerobic vs wild_production_aerobic, mutant_growth_aerobic vs mutant_production_aerobic

venn.diagram(x=list("wgaXwpa" = rownames(wild_growth_aerobic_X_wild_production_aerobic),
                    "mgaXmpa" = rownames(mutant_growth_aerobic_X_mutant_production_aerobic)), filename = "growth_X_production.png")

##wild_production_aerobic vs wild_production_anoxic, mutant_production_aerobic vs mutant_production_anoxic

venn.diagram(x=list("wpaXwpn" = rownames(wild_production_aerobic_X_wild_production_anoxic),
                    "mpaXmpn" = rownames(mutant_production_aerobic_X_mutant_production_anoxic)), filename = "aerobic_X_anoxic.png")
