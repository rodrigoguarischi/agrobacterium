# Agrobacterium

## Background

*Agrobacterium sp.* ATCC 31749 (also *Rhizobium sp.* ATCC 31749) is an alphaproteobacterium of the *Rhizobiaceae* family. 
This *Agrobacterium* is unique in its ability to synthesize, upon nitrogen exhaustion, a linear -1,3-glucan exopolysaccharide(EPS) known as curdlan.

## Purpose

First, we want to find the genes relevant to curdlan synthesis, nitrogen regulation and oxygen regulation, and we also want to know whether the genes relevant with oxygen have something with curdlan synthesis. 
Second, we want to research the function of *fnrN* gene in *Agrobacterium sp.* ATCC 31749. 

## Method

Samples of cell growth phase, curdlan-producing phase (aerobic) and curdlan-producing phase (anoxic treated) of both *Agrobacterium sp.* ATCC 31749 wild strain and ΔfnrN strain were collectecd to extract mRNA. 
Each sample was treated in duplicate. Illumina Hiseq4000 was used to complete the research.

## Sample description

**wild_growth_aerobic:** wild type strain, cell growth phase in aerobic culture condition.
**wild_production_aerobic:** wild type strain, curdlan production phase in aerobic culture condition.
**wild_production_anoxic:** wild type strain, curdlan production phase in anoxic culture condition.
**mutant_growth_aerobic:** ΔfnrN strain, cell growth phase in aerobic culture condition.
**mutant_production_aerobic:** ΔfnrN strain, curdlan production phase in aerobic culture condition.
**mutant_production_anoxic:** ΔfnrN strain, curdlan production phase in anoxic culture condition.

---

## Mission

1) **Read alignment will be performed using ~~Bowtie2~~ BWA using *Agrobacterium sp.* ATCC 31749's genome retrieved from NCBI. (Completed)**
2) Based on alignment files, quantification will be performed using HTSeq-count using RefSeq's genome annotation. Then the gene annotation analysis should include GO annotation, KEGG pathway annotation and COG annotation. Each gene GO annotation will be classified to BP, CC or MF. Each gene COG function classification and KEGG classification are also needed.
3) Differential expression analysis will be performed according your methods, then All comparisons will be made and displayed on a heatmap clustering genes by it's expression profiles.
4) **Each gene from RefSeq's annotation will be assigned to a COG category. (Completed)**
5) An enrichment analysis will be performed to determine enriched COGs on all comparisons between groups. And the enrichment analysis of GO terms and KEGG patnways on all comparisons between groups are also needed.
6) PCA analysis map, Cluster analysis map and Scatter plot map, Venn diagram of Differentially expressed genes on all comparisons between groups

### Data comparison analysis according to the following groups:

wild_growth_aerobic vs mutant_growth_aerobic
wild_production_aerobic vs mutant_production_aerobic
wild_production_anoxic vs mutant_production_anoxic
wild_growth_aerobic vs wild_production_aerobic
wild_production_aerobic vs wild_production_anoxic
mutant_growth_aerobic vs mutant_production_aerobic
mutant_production_aerobic vs mutant_production_anoxic

and

wild_growth_aerobic vs mutant_growth_aerobic, wild_production_aerobic vs mutant_production_aerobic, wild_production_anoxic vs mutant_production_anoxic
wild_growth_aerobic vs wild_production_aerobic, mutant_growth_aerobic vs mutant_production_aerobic
wild_production_aerobic vs wild_production_anoxic, mutant_production_aerobic vs mutant_production_anoxic
