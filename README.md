# rnaseqviewer

Automated visualization of RNA-seq differential expression analysis

# DEMO
Omics-overview
![sample_overview](https://user-images.githubusercontent.com/77435195/126021992-bcf85ab9-37ef-4409-adf0-d6d807abca12.png)

DEG_overview
![sample_pairwiseEBseq_viewer](https://user-images.githubusercontent.com/77435195/126033622-d33c24b8-14cd-4cd6-bd03-e32b1cd6c80a.png)

# Installation
```
install.packages("devtools")
devtools::install_github("Kan-E/rnaseqviewer")
```
# Usage
```
Omics_overview(Count_matrix = "normalized count data.txt")

DEG_overview(Count_matrix = "normalized count data.txt", 
             DEG_result = "result data of EBseq (or DEseq2).txt",
             Type = "EBseq",      #one of "EBseq" or "DEseq2"
             Species = "species",
             fdr = 0.05, fc = 2,  #fdr and fold change threshold
             Cond_1 = 3,          #sample number of condition_1
             Cond_2 = 3)          #sample number of condition_2
             
multiDEG_overview(Normalized_count_matrix = "normalized count data.txt", 
                  EBseq_result = "result data of EBseq.txt",
                  EBseq_condmeans = "result data of EBseq.txt",
                  Species = "species",               #human or mouse
                  fdr = 0.05, fc = 2, basemeam = 0,d  #fdr ,fold change, and basemean threshold
                  Cond_1 = 3,                        #sample number of condition_1
                  Cond_2 = 3,                        #sample number of condition_2
                  Cond_3 = 3)                        #sample number of condition_3
             , , ,
                           Species, Cond_1 = 3, Cond_2 = 3, Cond_3 = 3,
                           fdr = 0.05, fc = 2, basemean = 0
kmeansClustering(Count_matrix = "normalized count data.txt", 
                 Species = "species",
                 km,                  #number of k-means clustering
                 km_repeats,        #number of k-means runs to get a consensus k-means clustering
                 basemean_cutoff =0 )

GeneSetConversion(Gene_set = "gene list.txt")

AutoExtraction(Count_matrix = "normalized count data.txt", 
               Gene_set = "gene list.txt")
```
 
# Author
 
Kan Etoh
<kaneto@kumamoto-u.ac.jp>
