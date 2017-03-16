#VennDiagram (version 1.6.17)
library(VennDiagram)

#Carregando arquivos do DESeq2
mutant_growth_aerobicXmutant_production_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/mutant_growth_aerobicXmutant_production_aerobic.txt")
mutant_production_aerobicXmutant_production_anoxic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/mutant_production_aerobicXmutant_production_anoxic.txt")
wild_growth_aerobicXmutant_growth_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/wild_growth_aerobicXmutant_growth_aerobic.txt")
wild_growth_aerobicXwild_production_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/wild_growth_aerobicXwild_production_aerobic.txt")
wild_production_aerobicXmutant_production_aerobic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/wild_production_aerobicXmutant_production_aerobic.txt")
wild_production_aerobicXwild_production_anoxic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/wild_production_aerobicXwild_production_anoxic.txt")
wild_production_anoxicXmutant_production_anoxic <- read.delim("/mnt/work1/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/wild_production_anoxicXmutant_production_anoxic.txt")

#Gerando os diagramas de Venn

setwd("/work/agrobacterium/jhonatas/quant_gene_expression/htseq/DEG/deseq2/lfc_threshold/plots/")

##wild_growth_aerobic vs mutant_growth_aerobic, wild_production_aerobic vs mutant_production_aerobic, 
    ## wild_production_anoxic vs mutant_production_anoxic.

venn.diagram(x=list("Y1 vs F1" = rownames(wild_growth_aerobicXmutant_growth_aerobic),
                    "Y2 vs F2" = rownames(wild_production_aerobicXmutant_production_aerobic),
                    "Y3 vs F3" = rownames(wild_production_anoxicXmutant_production_anoxic)), filename = "wildXmutant.tiff", imagetype = "tiff", fill = c("green","blue","red"), col = "transparent")

##wild_growth_aerobic vs wild_production_aerobic, mutant_growth_aerobic vs mutant_production_aerobic

venn.diagram(x=list("Y1 vs Y2" = rownames(wild_growth_aerobicXwild_production_aerobic),
                    "F1 vs F2" = rownames(mutant_growth_aerobicXmutant_production_aerobic)), filename = "growthXproduction.tiff", imagetype = "tiff", fill = c("green","blue"), col = "transparent")

##wild_production_aerobic vs wild_production_anoxic, mutant_production_aerobic vs mutant_production_anoxic

venn.diagram(x=list("Y2 vs Y3" = rownames(wild_production_aerobicXwild_production_anoxic),
                    "F2 vs F3" = rownames(mutant_production_aerobicXmutant_production_anoxic)), filename = "aerobicXanoxic.tiff", imagetype = "tiff", fill = c("green","blue"), col = "transparent")
