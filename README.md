# rnaseqviewer

'rnaseqviewer' is an R package for automated visualization of RNA-seq differential expression analysis. It has three main functions, 'Omics_overview' for cluster analysis, and 'DEG_overview' and 'multiDEG_overview' for systematic pairwise and three-group comparison analysis of differentialy expressed genes, respectively.  It also has three useful functions, 'kmeansClustering' for k-means clustering analysis, 'AutoExtraction' for boxplotting gene sets, and 'vennd' for venn diagram analysis.  

# DEMO
Omics-overview
![sample_overview](https://user-images.githubusercontent.com/77435195/126021992-bcf85ab9-37ef-4409-adf0-d6d807abca12.png)

DEG_overview
![sample_pairwiseEBseq_viewer](https://user-images.githubusercontent.com/77435195/126033622-d33c24b8-14cd-4cd6-bd03-e32b1cd6c80a.png)

# Installation
```
install.packages("devtools")
devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
devtools::install_github("Kan-E/rnaseqviewer")
```
# Usage
```
Omics_overview(Count_matrix = "normalized count data.txt")

DEG_overview(Count_matrix,                            #normalized count data.txt 
             DEG_result,                              #result data of EBseq (or DEseq2).txt
             Type = "EBseq",                          #one of "EBseq" or "DEseq2"
             Species = NULL,                          #human or mouse (for enrichment analysis)
             fdr = 0.05, fc = 2, basemean = 0,        #fdr, fold change, and basemean threshold
             Cond_1 = 3,                              #sample number of condition_1
             Cond_2 = 3)                              #sample number of condition_2
             
multiDEG_overview(Normalized_count_matrix,            #normalized count data.txt 
                  EBseq_result,                       #result data of EBseq.txt
                  EBseq_condmeans,                    #result data of EBseq.txt"
                  Species = NULL,                     #human or mouse (for enrichment analysis)
                  fdr = 0.05, fc = 2, basemeam = 0,   #fdr ,fold change, and basemean threshold
                  Cond_1 = 3,                         #sample number of condition_1
                  Cond_2 = 3,                         #sample number of condition_2
                  Cond_3 = 3)                         #sample number of condition_3
             
kmeansClustering(Count_matrix,        #normalized count data.txt 
                 Species = NULL,      #Species for enrichment analysis
                 km,                  #number of k-means clustering
                 km_repeats = 10000,  #number of k-means runs to get a consensus k-means clustering
                 basemean =0 )        #basemean threshold

AutoExtraction(Count_matrix,        #normalized count data.txt 
               Gene_set_dir)        #directory including gene set txt files
               
vennd(gene_list_dir)

GeneSetConversion(Gene_set_dir)       #directory including gene set txt files

```
 
# Author
 
Kan Etoh
<kaneto@kumamoto-u.ac.jp>
