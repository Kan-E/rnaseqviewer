# rnaseqviewer

`rnaseqviewer` is an R package for automated visualization of RNA-seq differential expression analysis. It has functions for DEG analysis of single dataset and intergrated DEG analysis of multiple DEG analysis. `deseq2` and `ebseq` are functions for the DEG detection.  `Omics_overview`, `DEG_overview`, and `multiDEG_overview` are functions for systematic visualization of clustering and DEG analysis.  `vennd`, `ORA`, `AutoExtraction` and `int_heatmap` are functions for intergrated functional analysis.  It also has two useful functions, `kmeansClustering` for k-means clustering analysis, and `GeneSetConversion` for gene symbol conversion from human to mouse.  

![workflow](https://user-images.githubusercontent.com/77435195/139804204-f83030bc-8d69-4d74-a502-d990d94f6dec.png)

# DEMO
Omics-overview
![omics_overview](https://user-images.githubusercontent.com/77435195/132705815-11c55596-af12-439b-96cf-a961f39af2cf.png)

DEG_overview
![DEG_overview](https://user-images.githubusercontent.com/77435195/132705579-20bed45b-e9ce-4906-9e78-aaacea72d81a.png)

multiDEG_overview
![multiDEG_overview](https://user-images.githubusercontent.com/77435195/132705265-a87cb70c-cb8e-4d7e-bdc4-1f88c011cd3b.png)

# Installation
```
install.packages("devtools")
devtools::install_github("Kan-E/rnaseqviewer")
```
# Usage
```
Omics_overview(Count_matrix)                          #normalized count data.txt

DEG_overview(Count_matrix,                            #normalized count data.txt
             DEG_result,                              #result data of EBseq (or DEseq2).txt
             Species = NULL,                          #human or mouse (for enrichment analysis)
             fdr = 0.05, fc = 2, basemean = 0)        #fdr ,fold change, and basemean threshold

multiDEG_overview(Normalized_count_matrix,            #normalized count data.txt
                  EBseq_result,                       #result data of EBseq.txt
                  EBseq_condmeans,                    #result data of EBseq.txt"
                  Species = NULL,                     #human or mouse (for enrichment analysis)
                  fdr = 0.05, fc = 2, basemeam = 0)   #fdr ,fold change, and basemean threshold

kmeansClustering(Count_matrix,        #normalized count data.txt
                 Species = NULL,      #Species for enrichment analysis
                 km,                  #number of k-means clustering
                 km_repeats = 10000,  #number of k-means runs to get a consensus k-means clustering
                 basemean =0 )        #basemean threshold

AutoExtraction(Count_matrix,          #normalized count data.txt
               Gene_set_dir)          #directory including gene set txt files
               
int_heatmap(Count_matrix_dir,         #Directory including count matrix txt files
            Gene_set,                 #gene set txt file
            pre_zscoring = T)         #option for zscoring before integration of data sets

vennd(gene_list_dir)                  #directory including gene list txt files (up to 7 files)

GeneSetConversion(Gene_set_dir)       #directory including gene set txt files

deseq2(Row_count_matrix)              #Row count data.txt (NOT normalized count data)

ebseq(Row_count_matrix)              #Row count data.txt (NOT normalized count data)

```

# Reference
EBSeq (for ebseq)
- Ning Leng and Christina Kendziorski (2020). EBSeq: An R package for gene and isoform
  differential expression analysis of RNA-seq data. R package version 1.30.0.
  
DESeq2 (for deseq2)
- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

ggdendro (for dendrograms)
- Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro

clusterProfiler and DOSE (for enrichment analysis)
- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609

AnnotationDbi, org.Hs.eg.db and org.Mm.eg.db (for genome wide annotation)
- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi
- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.
- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.

genefilter (for z-score normalization)
- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.

ComplexHeatmap (for heatmap and k-means clustering)
- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.

ggplot2 and ggpubr (for boxplot and scater plot)
- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr

venn (for venn diagram analysis)
- Adrian Dusa (2021). venn: Draw Venn Diagrams. R package version 1.10. https://CRAN.R-project.org/package=venn

BioMart (for conversion of human gene symbol to mouse gene symbol)
- Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).

dplyr and tidyr (for data manipulation)
- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr

# Author

Kan Etoh
<kaneto@kumamoto-u.ac.jp>
