# Agrobacterium


## Background


*Agrobacterium sp.* ATCC 31749 (also *Rhizobium sp.* ATCC 31749) is an alphaproteobacterium of the *Rhizobiaceae* family.

This *Agrobacterium* is unique in its ability to synthesize, upon nitrogen exhaustion, a linear -1,3-glucan exopolysaccharide(EPS) known as curdlan.


## Purpose


First, we want to find the genes relevant to curdlan synthesis, nitrogen regulation and oxygen regulation, and we also want to know whether the genes relevant with oxygen have something with curdlan synthesis. 

Second, we want to research the function of *fnrN* gene in *Agrobacterium sp.* ATCC 31749. 


## Method


Samples of cell growth phase, curdlan-producing phase (aerobic) and curdlan-producing phase (anoxic treated) of both *Agrobacterium sp.* ATCC 31749 wild strain and ΔfnrN strain were collected to extract mRNA. 

Each sample was treated in duplicate. Illumina Hiseq4000 was used to complete the research.


## Sample description


**wild_growth_aerobic (Y1):** wild type strain, cell growth phase in aerobic culture condition.

**wild_production_aerobic (Y2):** wild type strain, curdlan production phase in aerobic culture condition.

**wild_production_anoxic (Y3):** wild type strain, curdlan production phase in anoxic culture condition.

**mutant_growth_aerobic (F1):** ΔfnrN strain, cell growth phase in aerobic culture condition.

**mutant_production_aerobic (F2):** ΔfnrN strain, curdlan production phase in aerobic culture condition.

**mutant_production_anoxic (F3):** ΔfnrN strain, curdlan production phase in anoxic culture condition.

---


## Mission


**1) Read alignment will be performed using Bowtie2 *Agrobacterium sp.* ATCC 31749's genome retrieved from NCBI. (Completed)**

**2) Based on alignment files, quantification will be performed using HTSeq-count using RefSeq's genome annotation. The quantification also will be performed using Salmon. (Completed)**

**3) Each gene from RefSeq's annotation will be assigned to a COG category. (Completed)**

**4) Each gene from RefSeq's annotation will be assigned to a GO term and each GO annotation will be classified to BP, CC or MF. (Completed)**

**5) Each gene from RefSeq's annotation will be assigned to a KEGG term and a KEGG pathway. (Completed)**

**6) Differential expression analysis will be performed with DESeq2 based on the comparison groups. (Completed)**

**7) An enrichment analysis will be performed to determine enriched COGs on all comparisons between groups. And the enrichment analysis of GO terms and KEGG pathways on all comparisons between groups are also needed. (Completed)**

**8) PCA analysis map and cluster analysis map (heatmap), Venn diagram of Differentially expressed genes on all comparisons between groups. (Completed)**


### Data comparison analysis according to the following groups:


wild_growth_aerobic (Y1) **vs** mutant_growth_aerobic (F1)

wild_production_aerobic (Y2) **vs** mutant_production_aerobic (F2)

wild_production_anoxic (Y3) **vs** mutant_production_anoxic (F3)

wild_growth_aerobic (Y1) **vs** wild_production_aerobic (Y2)

wild_production_aerobic (Y2) **vs** wild_production_anoxic (Y3)

mutant_growth_aerobic (F1) **vs** mutant_production_aerobic (F2)

mutant_production_aerobic (F2) **vs** mutant_production_anoxic (F3)


### Venn Diagrams of the following groups


wild_growth_aerobic (Y1) **vs** mutant_growth_aerobic (F1), wild_production_aerobic (Y2) **vs** mutant_production_aerobic (F2), wild_production_anoxic (Y3) **vs** mutant_production_anoxic (F3)

wild_growth_aerobic (Y1) **vs** wild_production_aerobic (Y2), mutant_growth_aerobic (F1) **vs** mutant_production_aerobic (F2)

wild_production_aerobic (Y2) **vs** wild_production_anoxic (Y3), mutant_production_aerobic (F2) **vs** mutant_production_anoxic (F3)
