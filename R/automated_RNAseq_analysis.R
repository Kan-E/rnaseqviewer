#' Automated extraction of gene set from count matrix and plot by boxplot
#'
#' @importFrom rstatix group_by
#' @importFrom rstatix add_xy_position
#' @importFrom rstatix tukey_hsd
#' @importFrom rstatix add_significance
#' @importFrom ggpubr ggboxplot
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 element_text
#' @importFrom tidyr gather
#' @importFrom dplyr %>%
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @examples library(rnaseqviewer)
#'
#' data(DEG_list)
#' dir.create("DEG_list")
#' write.table(DEG_list[1], file = "DEG_list/dataset1.txt", sep = "\t", quote = FALSE)
#' write.table(DEG_list[2], file = "DEG_list/dataset2.txt", sep = "\t", quote = FALSE)
#' vennd("DEG_list")
#'
#' data(Row_count_3conditions)
#' write.table(Row_count_3conditions, file = "Row_count_3conditions.txt", sep = "\t", quote = FALSE)
#'
#' AutoExtraction(Count_matrix = "Row_count_3conditions.txt", Gene_set_dir = "DEG_list/group_lists")
#'
#' @references Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
#' @references Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr
#' @references R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#' @references H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' @references Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr
#' @param Count_matrix count matrix txt file
#' @param Gene_set_dir Directory including gene set txt files
#' @export
#'
AutoExtraction <- function(Count_matrix, Gene_set_dir) {
  if(length(grep("/", Count_matrix) == 1)){
    dir_name <- paste0(gsub("/[^/]+$", "", Count_matrix), "/")
    file_name <- gsub(gsub("/[^/]+$", "", Count_matrix), "", Count_matrix)
  }else{
    dir_name <- ""
    file_name <- paste0("/", Count_matrix)
  }
  dir_name <- paste0(dir_name, "AutoExtraction")
  file_name <- gsub("\\..+$", "", file_name)
  dir.create(dir_name, showWarnings = F)
  dir_name_1 <- paste0(dir_name, paste0(file_name, "_boxplot"))
  dir_name_2 <- paste0(dir_name, paste0(file_name, "_test"))
  dir_name_3 <- paste0(dir_name, paste0(file_name, "_table"))
  dir_name_4 <- paste0(dir_name, paste0(file_name, "_heatmap"))
  dir.create(dir_name_1, showWarnings = F)
  dir.create(dir_name_2, showWarnings = F)
  dir.create(dir_name_3, showWarnings = F)
  dir.create(dir_name_4, showWarnings = F)

  All_data <- read.table(Count_matrix, header = T, row.names = 1, sep = "\t")
  group_files_full <- list.files(path = Gene_set_dir,
                               pattern = "*.txt", full.names = T)
  group_files <- list.files(path = Gene_set_dir,
                            pattern = "*.txt")
  group_dir <- gsub(group_files[1], "",group_files_full[1])
  group_files_full <- gsub(".txt", "", group_files_full)
  for (name in group_files_full) {
    data.file <- paste0(name, ".txt")
    print(data.file)
    group <- read.table(data.file, header = T, row.names = 1, sep = "\t")
    data <- merge(All_data, group, by = 0)
    data <- data[,1:(1 + ncol(All_data))]
    group.file <- paste(name, ".txt", sep = "")
    group.file <- gsub(group_dir, "", group.file)
    group.file <- paste(paste(dir_name_3, "/", sep = ""), group.file, sep = "")
    write.table(data, file = group.file, row.names = F,
                col.names = T, quote = F, sep = "\t")
    data <- read.table(group.file, header = T, sep = "\t")
    collist <- gsub("\\_.+$", "", colnames(data))
    collist <- unique(collist[-1])
    rowlist <- gsub("\\_.+$", "", data[,1])
    rowlist <- unique(rowlist)
    value <- NULL
    Row.names <- NULL
    data <- data %>% tidyr::gather(key = sample,
                                   value = value, -Row.names)
    data$sample <- gsub("\\_.+$", "", data$sample)
    data$Row.names <- as.factor(data$Row.names)
    data$sample <- factor(data$sample,levels=collist,ordered=TRUE)
    data$value <- as.numeric(data$value)
    stat.test <- data %>%
      group_by(Row.names) %>%
      tukey_hsd(value ~ sample) %>%
      add_significance("p.adj")
    stat.test <- stat.test %>% add_xy_position()
    stat.test
    if ((length(rowlist) > 81) && (length(rowlist) <= 200))
    {pdf_hsize <- 15
    pdf_wsize <- 15}
    if ((length(rowlist) > 64) && (length(rowlist) <= 81))
    {pdf_hsize <- 13.5
    pdf_wsize <- 13.5}
    if ((length(rowlist) > 49) && (length(rowlist) <= 64))
    {pdf_hsize <- 12
    pdf_wsize <- 12}
    if ((length(rowlist) > 36) && (length(rowlist) <= 49))
    {pdf_hsize <- 10.5
    pdf_wsize <- 10.5}
    if ((length(rowlist) > 25) && (length(rowlist) <= 36))
    {pdf_hsize <- 9
    pdf_wsize <- 9}
    if ((length(rowlist) > 16) && (length(rowlist) <= 25))
    {pdf_hsize <- 7.5
    pdf_wsize <- 7.5}
    if ((length(rowlist) > 12) && (length(rowlist) <= 16))
    {pdf_hsize <- 6
    pdf_wsize <- 6}
    if ((length(rowlist) > 9) && (length(rowlist) <= 12))
    {pdf_hsize <- 5
    pdf_wsize <- 6}
    if ((length(rowlist) > 6) && (length(rowlist) <= 9))
    {pdf_hsize <- 5
    pdf_wsize <- 4.5}
    if ((length(rowlist) > 4) && (length(rowlist) <= 6))
    {pdf_hsize <- 4
    pdf_wsize <- 6}
    if (length(rowlist) == 4)
    {pdf_hsize <- 4
    pdf_wsize <- 4}
    if (length(rowlist) == 3)
    {pdf_hsize <- 2
    pdf_wsize <- 6}
    if (length(rowlist) == 2)
    {pdf_hsize <- 2
    pdf_wsize <- 4}
    if (length(rowlist) == 1)
    {pdf_hsize <- 2
    pdf_wsize <- 2}
    if (length(rowlist) > 200)
    {pdf_hsize <- 30
    pdf_wsize <- 30}
    image.file2 <- paste(name, ".pdf", sep = "")
    image.file2 <- gsub(group_dir, "", image.file2)
    image.file2 <- paste(paste(dir_name_1, "/", sep = ""), image.file2, sep = "")
    pdf(image.file2, width = pdf_wsize, height = pdf_hsize)
    plot(ggboxplot(data, x = "sample", y = "value", fill = "sample",
                   facet.by = "Row.names", scales = "free", add = "jitter")
                  + theme(axis.text.x = element_text(size = 5),
                          axis.text.y = element_text(size = 10))
                  + scale_y_continuous(limits = c(0, NA)))
    dev.off()
    test.file <- paste(name, "_tukeyHSD.csv", sep = "")
    test.file <- gsub(group_dir, "", test.file)
    test.file <- paste(paste(dir_name_2, "/", sep = ""), test.file, sep = "")
    write.table(stat.test[, 1:10], file = test.file, sep = ",", quote = F, row.names = T)

    data<-read.table(group.file, header = T, row.names = 1, sep = "\t")
    data.z <- genescale(data, axis=1, method="Z")
    data.z <- na.omit(data.z)
    if(length(rowlist) <= 50){
    ht <- Heatmap(data.z, name = "z-score",
                  column_order = colnames(data.z),
                  clustering_method_columns = 'ward.D2',
                  show_row_names = T, show_row_dend = T)
    } else{
      ht <- Heatmap(data.z, name = "z-score",
                    column_order = colnames(data.z),
                    clustering_method_columns = 'ward.D2',
                    show_row_names = F, show_row_dend = T)
    }
    heatmap.file <- paste(name, '.pdf', sep = '')
    heatmap.file <- paste(paste(dir_name_4, "/", sep = ""), heatmap.file, sep = "")
    heatmap.file <- gsub(group_dir, "", heatmap.file)
    pdf(heatmap.file,width = 7,height = 10)
    print(ht)
    dev.off()
  }
}

#' Dendrogram, PCA, and mds
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggdendro dendro_data
#' @importFrom ggdendro segment
#' @importFrom ggdendro label
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ggplotify as.grob
#' @importFrom gridExtra grid.arrange
#' @importFrom gridExtra arrangeGrob
#' @importFrom utils read.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats cmdscale
#' @importFrom stats cor
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats na.omit
#' @importFrom stats prcomp
#' @references H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' @references Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro
#' @examples library(rnaseqviewer)
#' data(Row_count_3conditions)
#' write.table(Row_count_3conditions, file = "Row_count_3conditions.txt", sep = "\t", quote = FALSE)
#' ebseq("Row_count_3conditions.txt")
#' Omics_overview("Normalized_count_matrix_from_Row_count_3conditions_Cond1-vs-Cond2-vs-Cond3_EBseq.txt")
#' @param Count_matrix count matrix txt file
#' @param heatmap In the case of False, heatmap not shown.
#' @export
#'
Omics_overview <- function(Count_matrix, heatmap = TRUE){
  data <- read.table(Count_matrix, header = T, row.names = 1, sep = "\t")
  if(length(grep("/", Count_matrix) == 1)){
  dir_name <- paste0(gsub("/[^/]+$", "", Count_matrix), "/")
  file_name <- gsub(gsub("/[^/]+$", "", Count_matrix), "", Count_matrix)
  }else{
    dir_name <- ""
    file_name <- paste0("/", Count_matrix)
    }
  dir_name <- paste0(dir_name, "Omics_overview")
  dir.create(dir_name, showWarnings = F)
  file_name <- gsub("\\..+$", "", file_name)

  pca <- prcomp(data, scale. = T)
  label<- colnames(data)
  label<- gsub("\\_.+$", "", label)
  lab_x <- paste(summary(pca)$importance[2,1]*100,
                 "% of variance)", sep = "")
  lab_x <- paste("PC1 (", lab_x, sep = "")
  lab_y <- paste(summary(pca)$importance[2,2]*100,
                 "% of variance)", sep = "")
  lab_y <- paste("PC2 (", lab_y, sep = "")
  pca$rotation <- as.data.frame(pca$rotation)
  g1 <- ggplot(pca$rotation,aes(x=pca$rotation[,1],
                                y=pca$rotation[,2],
                                col=label, label = label)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA, size=0.5)) +
    xlab(lab_x) + ylab(lab_y) + geom_text_repel()  +
    theme(legend.position="none")
  pca_name <- paste0(dir_name, paste0(file_name, "_pca.pdf"))
  pdf(pca_name, height = 3, width = 3.5)
  print(g1)
  dev.off()

  rho <- cor(data,method="spearman")
  d <- dist(1-rho)
  mds <- as.data.frame(cmdscale(d))
  label<-colnames(data)
  label<-gsub("\\_.+$", "", label)
  g2 <- ggplot(mds, aes(x = mds[,1], y = mds[,2],
                        col = label, label = label)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA, size=0.5)) +
    xlab("dim 1") + ylab("dim 2") +
    geom_text_repel() + theme(legend.position="none")
  mds_name <- paste0(dir_name, paste0(file_name, "_mds.pdf"))
  pdf(mds_name, height = 3, width = 3.5)
  print(g2)
  dev.off()

  x <- NULL
  y <- NULL
  xend <- NULL
  yend <- NULL
  data.t <- t(data)
  hc <- hclust(dist(data.t), "ward.D2")
  dendr <- dendro_data(hc, type="rectangle")
  g3 <- ggplot() +
    geom_segment(data=segment(dendr),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=label(dendr),
              aes(x, y, label=label, hjust=0), size=2) +
    coord_flip()+ scale_y_reverse(expand=c(0.2, 0)) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  dg_name <- paste0(dir_name, paste0(file_name, "_wardD2.pdf"))
  pdf(dg_name, height = 3, width = 6)
  print(g3)
  dev.off()

  if(heatmap == TRUE){
  data.z <- genescale(data, axis=1, method="Z")
  data.z <- na.omit(data.z)
  ht <- as.grob(Heatmap(data.z, name = "z-score",
                        column_order = sort(colnames(data.z)),
                        clustering_method_columns = 'ward.D2',
                        show_row_names = F, show_row_dend = F))

  summary_name <- paste0(dir_name, paste0(file_name, "_summary.pdf"))
  pdf(summary_name, height = 9,width = 7)
  gridExtra::grid.arrange(gridExtra::arrangeGrob(g3, g1, g2, ncol = 1), ht, ncol = 2)
  dev.off()
  }
}


#' Visualization of pairwise DEG analysis
#'
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom ggpubr ggboxplot
#' @importFrom ggpubr ggmaplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 element_text
#' @importFrom tidyr gather
#' @importFrom dplyr %>%
#' @importFrom dplyr distinct
#' @importFrom AnnotationDbi select
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ggplotify as.grob
#' @importFrom gridExtra grid.arrange
#' @importFrom gridExtra arrangeGrob
#' @importFrom clusterProfiler compareCluster
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler enrichGO
#' @importFrom clusterProfiler gseGO
#' @importFrom clusterProfiler gseKEGG
#' @importFrom enrichplot dotplot
#' @importFrom enrichplot cnetplot
#' @importFrom enrichplot gseaplot2
#' @importFrom DOSE setReadable
#' @importFrom graphics barplot
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @import ggnewscale
#' @importFrom cowplot plot_grid
#' @references T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
#' @references Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609
#' @references Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0.
#' @references Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.
#' @references Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.
#' @references R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#' @references Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr
#' @references Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
#' @examples library(rnaseqviewer)
#'
#' data(Row_count_data)
#' write.table(Row_count_data, file = "Row_count_data.txt", sep = "\t", quote = FALSE)
#' deseq2("Row_count_data.txt")
#' DEG_overview(Count_matrix = "Normalized_count_matrix_from_Row_count_data_Cond1-vs-Cond2_DEseq2.txt",
#'              DEG_result = "result_of_Row_count_data_Cond1-vs-Cond2_DEseq2.txt",
#'              Species = "human", fc = 1.5)
#' @param Count_matrix Count matrix txt file
#' @param DEG_result result txt file of DEG analysis
#' @param Species Species
#' @param fdr Accepted false discovery rate for considering genes as differentially expressed
#' @param fc the fold change threshold. Only genes with a fold change >= fc and padj <= fdr are considered as significantly differentially expressed.
#' @param basemean basemean threshold
#' @export
#'
DEG_overview <- function(Count_matrix, DEG_result, Species = NULL,
                         fdr = 0.05, fc = 2, basemean = 0){
  data <- read.table(DEG_result,header = T, row.names = 1, sep = "\t")
  count <- read.table(Count_matrix,header=T, row.names = 1, sep = "\t")
  collist <- factor(gsub("\\_.+$", "", colnames(count)))
  vec <- c()
  for (i in 1:length(unique(collist))) {
    num <- length(collist[collist == unique(collist)[i]])
    vec <- c(vec, num)
  }
  Cond_1 <- vec[1]
  Cond_2 <- vec[2]
  Row.names <- NULL
  log2FoldChange <- NULL
  value <- NULL
  data <- merge(data,count, by=0)
  data <- dplyr::filter(data, apply(data[,8:(7 + Cond_1 + Cond_2)],1,mean) > basemean)
  dir_name <- gsub(".txt", "", Count_matrix)
  dir_name <- paste0(dir_name, paste0("_fc", fc))
  dir_name <- paste0(dir_name, paste0("_fdr", fdr))
  dir_name <- paste0(dir_name, paste0("_basemean", basemean))
  dir.create(dir_name, showWarnings = F)

  switch (colnames(data)[2],
          "PPEE" = Type <- "EBseq",
          "baseMean" = Type <- "DEseq2")
  if (!is.null(Species)){
  switch (Species,
          "mouse" = org <- org.Mm.eg.db,
          "human" = org <- org.Hs.eg.db)
  switch (Species,
          "mouse" = org_code <- "mmu",
          "human" = org_code <- "hsa")
  }
  if (Type == "EBseq"){
    data$padj <- data$PPEE
  }
  if (Type == "EBseq"){
    data$log2FoldChange <- -1 * log2(data$PostFC)
  }
  if (Type == "DEseq2"){
    data$log2FoldChange <- -1 * data$log2FoldChange
  }
  if (!is.null(Species)){
  my.symbols <- data$Row.names
  gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("ENTREZID", "SYMBOL", "UNIPROT"))
  colnames(gene_IDs) <- c("Row.names","ENTREZID","UNIPROT")
  data <- merge(data, gene_IDs, by="Row.names")
  data <- data %>% distinct(Row.names, .keep_all = T)
  }
  if (Type == "EBseq"){
  baseMean <- (data$C1Mean + data$C2Mean)*(1/2)
  data <- cbind(data, baseMean)
  }

  ##MA-plot
  m1 <- as.grob(ggmaplot(data, fdr = fdr, fc = fc, size = 0.4,
                         palette = c("#B31B21", "#1465AC", "darkgray"),
                         genenames = as.vector(data$Row.names),
                         legend = "top", top = 20,
                         font.label = c("bold", 6),font.legend = "bold",
                         font.main = "bold",
                         ggtheme = ggplot2::theme_minimal(),
                         select.top.method = "fc"))
  data2 <- dplyr::filter(data, abs(data$log2FoldChange) > log(fc, 2))
  if(nrow(data2) != 0){
  data2$group <- "upregulated"
  data2$group[data2$log2FoldChange < 0] <- "downregulated"
  data3 <- dplyr::filter(data2, abs(data2$padj) < fdr)

  ##heatmap
  data.z <- genescale(data3[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
  ht <- as.grob(Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                        clustering_method_columns = 'ward.D2',
                        show_row_names = F, show_row_dend = F))
  ma_name <- paste(paste(dir_name, "/", sep = ""),
                   "ma-heatmap.pdf", sep = "")
  pdf(ma_name, width = 7, height = 3.5)
  suppressWarnings(print(plot_grid(m1, ht, rel_widths = c(2, 1))))
  dev.off()

  if (!is.null(Species)){
  #dotplot
  universe <- AnnotationDbi::select(org,keys = rownames(count),
                                    keytype = "SYMBOL",
                                    columns = c("ENTREZID", "SYMBOL", "UNIPROT"))
  formula_res <- try(compareCluster(ENTREZID~group, data=data3,
                                    fun="enrichKEGG", organism=org_code,
                                    universe = universe), silent = T)
  if (class(formula_res) == "try-error") {
    formula_res <- NA
    p1 <- NULL
  } else{
  if ((length(as.data.frame(formula_res)) == 0) ||
      is.na(unique(as.data.frame(formula_res)$qvalue))) {
    p1 <- NULL
  } else{
    p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 7))
  }
  }
  ##cnetplot
  upgene <- data3[data3$log2FoldChange > log(fc, 2),]
  geneList_up <- upgene$log2FoldChange
  names(geneList_up) = as.character(upgene$ENTREZID)
  kk1 <- enrichKEGG(upgene$ENTREZID, organism =org_code,
                    pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  if(is.null(kk1)){
    cnet1 <- NULL
  } else cnet1 <- setReadable(kk1, org, 'ENTREZID')
  if (length(cnet1$ID) == 0) {
    p2 <- NULL
  } else{
  p2 <- as.grob(cnetplot(cnet1, foldChange=geneList_up,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE)+ guides(edge_color = "none"))
  keggenrich_name1 <- paste(paste(dir_name, "/", sep = ""),
                           "kegg_enrich_up.txt", sep = "")
  write.table(as.data.frame(cnet1), file = keggenrich_name1, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  downgene <- data3[data3$log2FoldChange < log(1/fc, 2),]
  geneList_down <- downgene$log2FoldChange
  names(geneList_down) = as.character(downgene$ENTREZID)
  kk2 <- enrichKEGG(downgene$ENTREZID, organism =org_code,
                    pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  if(is.null(kk2)){
    cnet2 <- NULL
  } else cnet2 <- setReadable(kk2, org, 'ENTREZID')
  if (length(cnet2$ID) == 0) {
    p3 <- NULL
  } else{
  p3 <- as.grob(cnetplot(cnet2, foldChange=geneList_down,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE)+ guides(edge_color = "none"))
  keggenrich_name2 <- paste(paste(dir_name, "/", sep = ""),
                           "kegg_enrich_down.txt", sep = "")
  write.table(as.data.frame(cnet2), file = keggenrich_name2, row.names = F, col.names = T, sep = "\t", quote = F)
  }

  ##GSEA plot
  data <- na.omit(data)
  geneList <- data$log2FoldChange
  names(geneList) = as.character(data$ENTREZID)
  geneList <- sort(geneList, decreasing = TRUE)
  kk3 <- gseKEGG(geneList = geneList, pvalueCutoff = 0.05,
                 organism = org_code, keyType = "kegg",
                 exponent = 1, eps = 0, pAdjustMethod = "none",
                 minGSSize = 50, maxGSSize = 500, by = "fgsea",
                 use_internal_data = FALSE, verbose = F)
  if (length(kk3$ID) == 0) {
    p4 <- NULL
  } else{
    kk3 <- setReadable(kk3, org, 'ENTREZID')
    if (length(kk3$ID) >= 5){
  p4 <- as.grob(gseaplot2(kk3, 1:5, pvalue_table = F))
    }else{
      p4 <- as.grob(gseaplot2(kk3, 1:length(kk3$ID), pvalue_table = F))
    }
  gsekegg_name <- paste(paste(dir_name, "/", sep = ""),
                        "gsekegg.txt", sep = "")
  write.table(as.data.frame(kk3), file = gsekegg_name, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  kegg_name <- paste(paste(dir_name, "/", sep = ""),
                     "kegg.pdf", sep = "")
  pdf(kegg_name, width = 11, height = 11)
  print(plot_grid(p1, p4, p2, p3, nrow = 2))
  dev.off()


  ##GOenrichment
  ##enrichment
  formula_res_go <- try(compareCluster(ENTREZID~group,
                                   data=data3,fun="enrichGO", OrgDb=org), silent =T)
  if (class(formula_res_go) == "try-error") {
    formula_res_go <- NA
    g1 <- NULL
  } else {
  if ((length(as.data.frame(formula_res_go)) == 0) ||
      is.na(unique(as.data.frame(formula_res_go)$qvalue))) {
    g1 <- NULL
  } else{
  g1 <- as.grob(dotplot(formula_res_go, color ="qvalue", font.size = 7))
  }
  }
  ##cnetplot
  go1 <- enrichGO(upgene$ENTREZID, OrgDb = org)
  if(is.null(go1)){
    cnet_go1 <- NULL
  } else cnet_go1 <- setReadable(go1, org, 'ENTREZID')
  if (length(cnet_go1$ID) == 0) {
    g2 <- NULL
  } else{
  g2 <- as.grob(cnetplot(cnet_go1, foldChange=geneList_up,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE)+ guides(edge_color = "none"))
  goenrich_name1 <- paste(paste(dir_name, "/", sep = ""),
                         "go_enrich_up.txt", sep = "")
  write.table(as.data.frame(cnet_go1), file = goenrich_name1,
              row.names = F, col.names = T, sep = "\t", quote = F)
  }
  go2 <- enrichGO(downgene$ENTREZID, OrgDb = org)
  if(is.null(go2)){
    cnet_go2 <- NULL
  } else cnet_go2 <- setReadable(go2, org, 'ENTREZID')
  if (length(cnet_go2$ID) == 0) {
    g3 <- NULL
  } else{
  g3 <- as.grob(cnetplot(cnet_go2, foldChange=geneList_down,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE)+ guides(edge_color = "none"))
  goenrich_name2 <- paste(paste(dir_name, "/", sep = ""),
                         "go_enrich_down.txt", sep = "")
  write.table(as.data.frame(cnet_go2), file = goenrich_name2,
              row.names = F, col.names = T, sep = "\t", quote = F)
  }
  ##GSEA plot
  data <- na.omit(data)
  geneList_2 <- data$log2FoldChange
  names(geneList_2) = as.character(data$ENTREZID)
  geneList_2 <- sort(geneList_2, decreasing = TRUE)
  go3 <- gseGO(geneList = geneList_2, pvalueCutoff = 0.05,
               OrgDb = org, exponent = 1, eps = 0,
               pAdjustMethod = "none", minGSSize = 50,
               maxGSSize = 500, by = "fgsea", verbose = F)
  if (length(go3$ID) == 0) {
    g4 <- NULL
  } else{
    go3 <- setReadable(go3, org, 'ENTREZID')
    if (length(go3$ID) >= 5) {
  g4 <- as.grob(gseaplot2(go3, 1:5, pvalue_table = F))
    }else{
      g4 <- as.grob(gseaplot2(go3, 1:length(go3$ID), pvalue_table = F))
      }
        gsego_name <- paste(paste(dir_name, "/", sep = ""),
                      "gseGO.txt", sep = "")
  write.table(as.data.frame(go3), file = gsego_name, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  go_name <- paste(paste(dir_name, "/", sep = ""),
                   "GO.pdf", sep = "")
  pdf(go_name, width = 11, height = 11)
  print(plot_grid(g1, g4, g2, g3, nrow = 2))
  dev.off()
  }

  #FC上位50下位50をboxplot
  data4 <- data3[sort(data3$log2FoldChange, decreasing = T, index=T)$ix,]
  up_all <- dplyr::filter(data4, log2FoldChange > 0)
  up50 <- up_all[1:50,8:(7 + Cond_1 + Cond_2)]
  collist <- gsub("\\_.+$", "", colnames(up50))
  collist <- unique(collist[-1])
  up50$Row.names <- up_all[1:50,]$Row.names
  up50 <- up50 %>% gather(key=sample, value=value,-Row.names)
  up50$sample <- gsub("\\_.+$", "", up50$sample)
  up50$Row.names <- as.factor(up50$Row.names)
  up50$sample <- factor(up50$sample,levels=collist,ordered=TRUE)
  up50$value <- as.numeric(up50$value)
  data4 <- data3[sort(data3$log2FoldChange, decreasing = F, index=T)$ix,]
  down_all <- dplyr::filter(data4, log2FoldChange < 0)
  down50 <- down_all[1:50,8:(7 + Cond_1 + Cond_2)]
  down50$Row.names <- down_all[1:50,]$Row.names
  down50 <- down50 %>% gather(key=sample, value=value,-Row.names)
  down50$sample <- gsub("\\_.+$", "", down50$sample)
  down50$Row.names <- as.factor(down50$Row.names)
  down50$sample <- factor(down50$sample,levels=collist,ordered=TRUE)
  down50$value <- as.numeric(down50$value)

  data5 <- data4[,8:(7 + Cond_1 + Cond_2)]
  rownames(data5) <- data4$Row.names
  deg_name <- paste(paste(dir_name, "/", sep = ""),
                    "DEG_count.txt", sep = "")
  write.table(data5, file = deg_name, row.names = T, col.names = T, sep = "\t", quote = F)
  rownames(up_all) <- up_all$Row.names
  degup_name <- paste(paste(dir_name, "/", sep = ""),
                    "DEG_count_up.txt", sep = "")
  write.table(up_all[,8:(7 + Cond_1 + Cond_2)], file = degup_name, row.names = T, col.names = T, sep = "\t", quote = F)
  rownames(down_all) <- down_all$Row.names
  degdown_name <- paste(paste(dir_name, "/", sep = ""),
                      "DEG_count_down.txt", sep = "")
  write.table(down_all[,8:(7 + Cond_1 + Cond_2)], file = degdown_name, row.names = T, col.names = T, sep = "\t", quote = F)


  degtop_name <- paste(paste(dir_name, "/", sep = ""),
                       "top_DEG.pdf", sep = "")
  pdf(degtop_name,height = 10, width = 10)
  suppressWarnings(plot(ggpubr::ggboxplot(up50, x = "sample", y = "value",
                         fill = "sample", facet.by = "Row.names",
                         scales = "free", add = "jitter",
                         xlab = "gene", ylab = "Normalized_count")+
         ggplot2::theme(axis.text.x= ggplot2::element_text(size = 5),
                        axis.text.y= ggplot2::element_text(size = 10)) +
         ggplot2::scale_y_continuous(limits = c(0, NA)) +
    scale_fill_manual(values=c("gray", "#ff8082"))))
  suppressWarnings(plot(ggpubr::ggboxplot(down50, x = "sample", y = "value",
                         fill = "sample", facet.by = "Row.names",
                         scales = "free", add = "jitter",
                         xlab = "gene", ylab = "Normalized_count")+
         ggplot2::theme(axis.text.x= ggplot2::element_text(size = 5),
                        axis.text.y= ggplot2::element_text(size = 10)) +
         ggplot2::scale_y_continuous(limits = c(0, NA)) +
         scale_fill_manual(values=c("gray", "#ff8082"))))
  dev.off()
  } else{
    ma_name <- paste(paste(dir_name, "/", sep = ""),
                     "maplot.pdf", sep = "")
    pdf(ma_name, width = 4.5, height = 3.5)
    suppressWarnings(print(plot_grid(m1)))
    dev.off()
  }
}


#' Kmeans clustering analysis
#'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 element_text
#' @importFrom tidyr gather
#' @importFrom dplyr %>%
#' @importFrom dplyr distinct
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap row_dend
#' @importFrom ComplexHeatmap row_order
#' @importFrom graphics barplot
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @import ggnewscale
#' @references R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#' @param Count_matrix Directory of count matrix txt file
#' @param km number of k-means clustering
#' @param km_repeats number of k-means runs to get a consensus k-means clustering
#' @param basemean basemean threshold
#' @export
#'
kmeansClustring <- function(Count_matrix, km, km_repeats=10000,
                            basemean = 0){
  dir_name <- gsub("\\..+$", "", Count_matrix)
  dir_name <- paste(dir_name, paste("_km", km, sep = ""), sep = "")
  dir.create(dir_name, showWarnings = F)
  RNAseq <- read.table(Count_matrix, header=T, row.names = 1, sep = "\t")
  RNAseq2 <- dplyr::filter(RNAseq, apply(RNAseq,1,mean) > basemean)
  data.z <- genescale(RNAseq2, axis = 1, method = "Z")
  ht <- Heatmap(data.z, name = "z-score",
                column_order = colnames(data.z),
              clustering_method_columns = 'ward.D2',
              row_km= km, cluster_row_slices = F, row_km_repeats = km_repeats,
              show_row_names = F)
  ht <- draw(ht)
  r.dend <- row_dend(ht)
  rcl.list <- row_order(ht)
  lapply(rcl.list, function(x) length(x))
  Cluster <- NULL
  for (i in 1:length(row_order(ht))){ if (i == 1) {
    clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")} else {
      clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)}}
  cluster.file <- paste0(dir_name, "/gene_cluster.txt")
  write.table(out, file= cluster.file, sep="\t", quote=F, row.names=FALSE)
  out2 <- read.table(cluster.file, header = T, row.names = 1, sep = "\t")
  clusterlist <- unique(out2$Cluster)
  for (cluster_name in clusterlist) {
  clusterCount <- dplyr::filter(out2, Cluster == cluster_name)
  clusterCount <- merge(clusterCount, RNAseq2, by=0)
  clusterCount <- clusterCount[,-2]
  table.file <- paste0(paste0(dir_name,"/Count_"), paste0(cluster_name, ".txt"))
  write.table(clusterCount, file= table.file, sep="\t", quote=F, row.names = F)
  }
  heatmap.file <- paste0(dir_name, "/heatmap.pdf")
  pdf(heatmap.file, width = 3, height = 4)
  print(ht)
  dev.off()
}



#' Gene symbol conversion from human to mouse
#'
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getLDS
#' @examples library(rnaseqviewer)
#'
#' data(DEG_list)
#' dir.create("DEG_list")
#' write.table(DEG_list[1], file = "DEG_list/dataset1.txt", sep = "\t", quote = FALSE)
#' write.table(DEG_list[2], file = "DEG_list/dataset2.txt", sep = "\t", quote = FALSE)
#' GeneSetConversion("DEG_list")
#' @references Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).
#' @param Gene_set_dir Directory including Gene_set txt files
#' @export
#'
GeneSetConversion <- function(Gene_set_dir) {
  dir_name <- Gene_set_dir
  dir_name_1 <- paste(dir_name, "_convert", sep = "")
  dir.create(dir_name_1, showWarnings = F)
  gene_set_files <- list.files(path = Gene_set_dir,
                               pattern = "*.txt")
  gene_set_files <- gsub(".txt", "", gene_set_files)
  gene_files_full <- list.files(path = Gene_set_dir,
                                pattern = "*.txt", full.names = T)
  gene_files_full <- gsub(".txt", "", gene_files_full)
  gene_set_dir <- gsub(gene_set_files[1], "", gene_files_full[1])
  gene_files_full <- gsub(".txt", "", gene_files_full)

  for (name in gene_files_full) {
    data.file <- paste(name, ".txt", sep = "")
    print(data.file)
    genes <- read.table(data.file, sep = "\t", row.names = 1)
    genes <- rownames(genes)
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://www.ensembl.org")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
    genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                   values = genes ,mart = mouse,
                   attributesL = c("hgnc_symbol"),
                   martL = human, uniqueRows=T)
    gene.file <- gsub(gene_set_dir, "", data.file)
    gene.file2 <- paste(paste(dir_name_1, "/", sep = ""), gene.file, sep = "")
    write.table(genes, file = gene.file2, sep="\t", quote=F, row.names = F)
  }
}

