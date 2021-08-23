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
#' @importFrom readr write_csv
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @param Count_matrix count matrix txt file
#' @param Gene_set Directory of gene set
#' @export
#'
AutoExtraction <- function(Count_matrix, Gene_set) {
  dir_name <- gsub("\\..+$", "", Count_matrix)
  dir_name_1 <- paste(dir_name, "_boxplot", sep = "")
  dir_name_2 <- paste(dir_name, "_boxplot_noAsterisk", sep = "")
  dir_name_3 <- paste(dir_name, "_test", sep = "")
  dir_name_4 <- paste(dir_name, "_table", sep = "")
  dir_name_5 <- paste(dir_name, "_heatmap", sep = "")
  dir.create(dir_name_1, showWarnings = F)
  dir.create(dir_name_2, showWarnings = F)
  dir.create(dir_name_3, showWarnings = F)
  dir.create(dir_name_4, showWarnings = F)
  dir.create(dir_name_5, showWarnings = F)
  All_data <- read.table(Count_matrix, header = T, row.names = 1)
  group_files_full <- list.files(path = Gene_set,
                               pattern = "*.txt", full.names = T)
  group_files <- list.files(path = Gene_set,
                            pattern = "*.txt")
  group_dir <- gsub(group_files[1], "",group_files_full[1])
  group_files_full <- gsub(".txt", "", group_files_full)
  for (name in group_files_full) {
    data.file <- paste(name, ".txt", sep = "")
    print(data.file)
    group <- read.table(data.file, header = F,
                           row.names = 1, skip = 2)
    data <- merge(group, All_data, by = 0)
    group.file <- paste(name, ".csv", sep = "")
    group.file <- gsub(group_dir, "", group.file)
    group.file <- paste(paste(dir_name_4, "/", sep = ""), group.file, sep = "")
    write.table(data, file = group.file, row.names = F,
                col.names = T, quote = F, sep = ",")
    data <- read.csv(group.file, header = T)
    collist <- gsub("\\_.+$", "", colnames(data))
    collist <- unique(collist[-1])
    rowlist <- gsub("\\_.+$", "", data[,1])
    rowlist <- unique(rowlist)
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
    if ((length(rowlist) > 81) && (length(rowlist) <= 100)) pdf_size <- 15
    if ((length(rowlist) > 64) && (length(rowlist) <= 81)) pdf_size <- 13.5
    if ((length(rowlist) > 49) && (length(rowlist) <= 64)) pdf_size <- 12
    if ((length(rowlist) > 36) && (length(rowlist) <= 49)) pdf_size <- 10.5
    if ((length(rowlist) > 25) && (length(rowlist) <= 36)) pdf_size <- 9
    if ((length(rowlist) > 16) && (length(rowlist) <= 25)) pdf_size <- 7.5
    if ((length(rowlist) > 12) && (length(rowlist) <= 16)) pdf_size <- 6
    if ((length(rowlist) > 9) && (length(rowlist) <= 12)) pdf_size <- 6
    if ((length(rowlist) > 6) && (length(rowlist) <= 9)) pdf_size <- 5
    if ((length(rowlist) > 2) && (length(rowlist) <= 6)) pdf_size <- 4
    if (length(rowlist) == 1) pdf_size <- 3
    if (length(rowlist) > 100) pdf_size <- 16.5
    image.file <- paste(name, ".pdf", sep = "")
    image.file <- gsub(group_dir, "", image.file)
    image.file <- paste(paste(dir_name_1, "/", sep = ""), image.file, sep = "")
    pdf(image.file, width = pdf_size, height = pdf_size)
    plot(ggboxplot(data, x = "sample", y = "value", fill = "sample",
                   facet.by = "Row.names", scales = "free", add = "jitter")
                    + stat_pvalue_manual(stat.test, hide.ns = T,
                                         tip.length = 0.01)
                    + theme(axis.text.x = element_text(size = 5),
                            axis.text.y = element_text(size = 10)))
    dev.off()
    image.file2 <- paste(name, "_noAsterisk.pdf", sep = "")
    image.file2 <- gsub(group_dir, "", image.file2)
    image.file2 <- paste(paste(dir_name_2, "/", sep = ""), image.file2, sep = "")
    pdf(image.file2, width = pdf_size, height = pdf_size)
    plot(ggboxplot(data, x = "sample", y = "value", fill = "sample",
                   facet.by = "Row.names", scales = "free", add = "jitter")
                  + theme(axis.text.x = element_text(size = 5),
                          axis.text.y = element_text(size = 10))
                  + scale_y_continuous(limits = c(0, NA)))
    dev.off()
    test.file <- paste(name, "_tukeyHSD.csv", sep = "")
    test.file <- gsub(group_dir, "", test.file)
    test.file <- paste(paste(dir_name_3, "/", sep = ""), test.file, sep = "")
    write_csv(stat.test[, 1:10], file = test.file)

    data<-read.csv(group.file, header = T, row.names = 1)
    data.z <- genescale(data, axis=1, method="Z")
    data.z <- na.omit(data.z)
    ht <- Heatmap(data.z, name = "z-score",
                  clustering_method_columns = 'ward.D2',
                  show_row_names = T, show_row_dend = T)
    heatmap.file <- paste(name, '.pdf', sep = '')
    heatmap.file <- paste(paste(dir_name_5, "/", sep = ""), heatmap.file, sep = "")
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
#' @param Count_matrix count matrix txt file
#' @export
#'
Omics_overview <- function(Count_matrix){
  data <- read.table(Count_matrix, header = T, row.names = 1)
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
  pca_name <- gsub("\\..+$", "", Count_matrix)
  pca_name <- paste(pca_name, "_pca.pdf", sep = "")
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
    xlab("Leading LogFC (dim 1)") + ylab("Leading LogFC (dim 2)") +
    geom_text_repel() + theme(legend.position="none")
  mds_name <- gsub("\\..+$", "", Count_matrix)
  mds_name <- paste(mds_name, "_mds.pdf", sep = "")
  pdf(mds_name, height = 3, width = 3.5)
  print(g2)
  dev.off()

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
  dg_name <- gsub("\\..+$", "", Count_matrix)
  dg_name <- paste(dg_name, "_wardD2.pdf", sep = "")
  pdf(dg_name, height = 3, width = 6)
  print(g3)
  dev.off()

  data.z <- genescale(data, axis=1, method="Z")
  data.z <- na.omit(data.z)
  ht <- as.grob(Heatmap(data.z, name = "z-score",
                        column_order = sort(colnames(data.z)),
                        clustering_method_columns = 'ward.D2',
                        show_row_names = F, show_row_dend = F))

  summary_name <- gsub("\\..+$", "", Count_matrix)
  summary_name <- paste(summary_name, "_summary.pdf", sep = "")
  pdf(summary_name, height = 9,width = 7)
  gridExtra::grid.arrange(gridExtra::arrangeGrob(g3, g1, g2, ncol = 1), ht, ncol = 2)
  dev.off()
}


#' Visualization of pairwise DEG analysis
#'
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom rstatix group_by
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
#' @param Count_matrix Count matrix txt file
#' @param DEG_result result txt file of DEG analysis
#' @param Type one of "EBseq" or "DEseq2"
#' @param Species Species
#' @param Cond_1 Sumple number of condition_1
#' @param Cond_2 Sumple number of condition_2
#' @param fdr Accepted false discovery rate for considering genes as differentially expressed
#' @param fc the fold change threshold. Only genes with a fold change >= fc and padj <= fdr are considered as significantly differentially expressed.
#' @export
#'
DEG_overview <- function(Count_matrix, DEG_result, Type = "EBseq",
                                 Species, Cond_1 = 3, Cond_2 = 3,
                                 fdr = 0.05, fc = 2){
  data <- read.table(DEG_result,header = T, row.names = 1)
  count <- read.table(Count_matrix,header=T, row.names = 1)
  data <- merge(data,count, by=0)

  dir_name <- gsub(".txt", "", Count_matrix)
  dir_name <- paste(dir_name, paste("_fc", fc, sep = ""), sep = "")
  dir_name <- paste(dir_name, paste("_fdr", fdr, sep = ""), sep = "")
  dir.create(dir_name, showWarnings = F)

  switch (Species,
          "mouse" = org <- org.Mm.eg.db,
          "human" = org <- org.Hs.eg.db)
  switch (Species,
          "mouse" = org_code <- "mmu",
          "human" = org_code <- "hsa")
  if (Type == "EBseq"){
    data$padj <- data$PPEE
  }
  if (Type == "EBseq"){
    data$log2FoldChange <- -1 * log2(data$PostFC)
  }
  my.symbols <- data$Row.names
  gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("ENTREZID", "SYMBOL", "UNIPROT"))
  colnames(gene_IDs) <- c("Row.names","ENTREZID","UNIPROT")
  data <- merge(data, gene_IDs, by="Row.names")
  data <- data %>% distinct(Row.names, .keep_all = T)
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
  data2 <- data[abs(data$log2FoldChange) > log(fc, 2),]
  if(nrow(data2) != 0){
  data2$group <- "upregulated"
  data2$group[data2$log2FoldChange < 0] <- "downregulated"
  data3 <- data2[abs(data2$padj) < fdr,]

  ##heatmap
  data.z <- genescale(data3[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
  ht <- as.grob(Heatmap(data.z, name = "z-score",
                        clustering_method_columns = 'ward.D2',
                        show_row_names = F, show_row_dend = T))
  ma_name <- paste(paste(dir_name, "/", sep = ""),
                   "ma-heatmap.pdf", sep = "")
  pdf(ma_name, width = 7, height = 3.5)
  print(plot_grid(m1, ht, rel_widths = c(2, 1)))
  dev.off()

  #dotplot
  universe <- AnnotationDbi::select(org,keys = rownames(count),
                                    keytype = "SYMBOL",
                                    columns = c("ENTREZID", "SYMBOL", "UNIPROT"))
  formula_res <- try(compareCluster(ENTREZID~group, data=data3,
                                    fun="enrichKEGG", organism=org_code,
                                    universe = universe), silent = T)
  if (class(formula_res) == "try-error") {
    formula_res <- NA
  }
  if ((length(as.data.frame(formula_res)) == 0) ||
      is.na(unique(as.data.frame(formula_res)$qvalue)) ||
      is.na(formula_res)) {
    p1 <- NULL
  } else{
    p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 9))
    keggenrich_name <- paste(paste(dir_name, "/", sep = ""),
          "kegg_enrich.txt", sep = "")
    write.table(as.data.frame(formula_res), file = keggenrich_name, row.names = F, col.names = T, sep = "\t", quote = F)
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
                         cex_category = 0.5, colorEdge = TRUE))
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
                         cex_category = 0.5, colorEdge = TRUE))
  }

  ##GSEA plot
  data_uniprot <- na.omit(data)
  data_uniprot <- data_uniprot %>% distinct(UNIPROT, .keep_all = T)
  geneList <- data_uniprot$log2FoldChange
  names(geneList) = as.character(data_uniprot$UNIPROT)
  geneList <- sort(geneList, decreasing = TRUE)
  kk2 <- gseKEGG(geneList = geneList, pvalueCutoff = 0.05,
                 organism = org_code, keyType = "uniprot",
                 exponent = 1, eps = 0, pAdjustMethod = "none",
                 minGSSize = 50, maxGSSize = 500, by = "fgsea",
                 use_internal_data = FALSE)
  if (length(kk2$ID) == 0) {
    p4 <- NULL
  } else{
  p4 <- as.grob(gseaplot2(kk2, 1:6, pvalue_table = F))
  gsekegg_name <- paste(paste(dir_name, "/", sep = ""),
                        "gsekegg.txt", sep = "")
  write.table(as.data.frame(kk2), file = gsekegg_name, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  kegg_name <- paste(paste(dir_name, "/", sep = ""),
                     "kegg.pdf", sep = "")
  pdf(kegg_name, width = 14, height = 14.0)
  print(plot_grid(p1, p4, p2, p3, nrow = 2))
  dev.off()


  ##GOenrichment
  ##enrichment
  formula_res_go <- try(compareCluster(ENTREZID~group,
                                   data=data3,fun="enrichGO", OrgDb=org), silent =T)
  if (class(formula_res_go) == "try-error") {
    formula_res_go <- NA
  }
  if ((length(as.data.frame(formula_res_go)) == 0) ||
      is.na(unique(as.data.frame(formula_res_go)$qvalue)) ||
      is.na(formula_res_go)) {
    g1 <- NULL
  } else{
  g1 <- as.grob(dotplot(formula_res_go, color ="qvalue", font.size = 9))
  goenrich_name <- paste(paste(dir_name, "/", sep = ""),
                         "go_enrich.txt", sep = "")
  write.table(as.data.frame(formula_res_go), file = goenrich_name,
              row.names = F, col.names = T, sep = "\t", quote = F)
  }
  ##cnetplot
  go1 <- enrichGO(upgene$ENTREZID, OrgDb = org,
                  pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  if(is.null(go1)){
    cnet_go1 <- NULL
  } else cnet_go1 <- setReadable(go1, org, 'ENTREZID')
  if (length(cnet_go1$ID) == 0) {
    g2 <- NULL
  } else{
  g2 <- as.grob(cnetplot(cnet_go1, foldChange=geneList_up,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE))
  }
  go2 <- enrichGO(downgene$ENTREZID, OrgDb = org,
                  pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  if(is.null(go2)){
    cnet_go2 <- NULL
  } else cnet_go2 <- setReadable(go2, org, 'ENTREZID')
  if (length(cnet_go2$ID) == 0) {
    g3 <- NULL
  } else{
  g3 <- as.grob(cnetplot(cnet_go2, foldChange=geneList_down,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE))
  }
  ##GSEA plot
  data <- na.omit(data)
  geneList <- data$log2FoldChange
  names(geneList) = as.character(data$ENTREZID)
  geneList <- sort(geneList, decreasing = TRUE)
  go3 <- gseGO(geneList = geneList, pvalueCutoff = 0.05,
               OrgDb = org, exponent = 1, eps = 0,
               pAdjustMethod = "none", minGSSize = 50,
               maxGSSize = 500, by = "fgsea")
  if (length(go3$ID) == 0) {
    g4 <- NULL
  } else{
  g4 <- as.grob(gseaplot2(kk2, 1:6, pvalue_table = F))
  gsego_name <- paste(paste(dir_name, "/", sep = ""),
                      "gseGO.txt", sep = "")
  write.table(as.data.frame(kk2), file = gsego_name, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  go_name <- paste(paste(dir_name, "/", sep = ""),
                   "GO.pdf", sep = "")
  pdf(go_name, width = 14, height = 14.0)
  print(plot_grid(g1, g4, g2, g3, nrow = 2))
  dev.off()

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

  degtop_name <- paste(paste(dir_name, "/", sep = ""),
                       "top_DEG.pdf", sep = "")
  pdf(degtop_name,height = 10, width = 10)
  plot(ggpubr::ggboxplot(up50, x = "sample", y = "value",
                         fill = "sample", facet.by = "Row.names",
                         scales = "free", add = "jitter",
                         xlab = "gene", ylab = "TPM")+
         ggplot2::theme(axis.text.x= ggplot2::element_text(size = 5),
                        axis.text.y= ggplot2::element_text(size = 10)) +
         ggplot2::scale_y_continuous(limits = c(0, NA)) +
    scale_fill_manual(values=c("gray", "#ff8082")))
  plot(ggpubr::ggboxplot(down50, x = "sample", y = "value",
                         fill = "sample", facet.by = "Row.names",
                         scales = "free", add = "jitter",
                         xlab = "gene", ylab = "TPM")+
         ggplot2::theme(axis.text.x= ggplot2::element_text(size = 5),
                        axis.text.y= ggplot2::element_text(size = 10)) +
         ggplot2::scale_y_continuous(limits = c(0, NA)) +
         scale_fill_manual(values=c("gray", "#ff8082")))
  dev.off()
  } else{
    ma_name <- paste(paste(dir_name, "/", sep = ""),
                     "maplot.pdf", sep = "")
    pdf(ma_name, width = 4.5, height = 3.5)
    print(plot_grid(m1))
    dev.off()
  }
}


#' Kmeans clustering analysis
#'
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom rstatix group_by
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
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap row_dend
#' @importFrom ComplexHeatmap row_order
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
#' @param Count_matrix Directory of count matrix txt file
#' @param Species Species
#' @param km number of k-means clustering
#' @param km_repeats number of k-means runs to get a consensus k-means clustering
#' @param basemean_cutoff basemean cutoff
#' @export
#'
kmeansClustring <- function(Count_matrix, Species, km, km_repeats,
                            basemean_cutoff = 0){
  dir_name <- gsub("\\..+$", "", Count_matrix)
  dir_name <- paste(dir_name, paste("_km", km, sep = ""), sep = "")
  dir.create(dir_name, showWarnings = F)
  RNAseq <- read.table(Count_matrix, header=T, row.names = 1)
  RNAseq2 <- dplyr::filter(RNAseq, apply(RNAseq,1,mean) > basemean_cutoff)
  data.z <- genescale(RNAseq2, axis = 1, method = "Z")
  ht <- Heatmap(data.z, name = "z-score",
              clustering_method_columns = 'ward.D2',
              row_km= km, cluster_row_slices = F, row_km_repeats = km_repeats,
              show_row_names = F)
  ht <- draw(ht)
  r.dend <- row_dend(ht)
  rcl.list <- row_order(ht)
  lapply(rcl.list, function(x) length(x))
  for (i in 1:length(row_order(ht))){ if (i == 1) {
    clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")} else {
      clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)}}
  cluster.file <- paste0(dir_name, "/gene_cluster.txt")
  write.table(out, file= cluster.file, sep="\t", quote=F, row.names=FALSE)
  out2 <- read.table(cluster.file, header = T, row.names = 1)
  clusterlist <- unique(out2$Cluster)
  for (cluster_name in clusterlist) {
  clusterTPM <- dplyr::filter(out2, Cluster == cluster_name)
  clusterTPM <- merge(clusterTPM, RNAseq2, by=0)
  clusterTPM <- clusterTPM[,-2]
  table.file <- paste0(paste0(dir_name,"/TPM_"), paste0(cluster_name, ".txt"))
  write.table(clusterTPM, file= table.file, sep="\t", quote=F, row.names = F)
  }

  switch (Species,
          "mouse" = org <- org.Mm.eg.db,
          "human" = org <- org.Hs.eg.db)
  switch (Species,
          "mouse" = org_code <- "mmu",
          "human" = org_code <- "hsa")
  my.symbols <- out
  gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                    keytype = "SYMBOL",
                                    columns = c("ENTREZID", "SYMBOL"))
  colnames(gene_IDs) <- c("GeneID","ENTREZID")
  data2 <- merge(out, gene_IDs, by="GeneID")
  formula_res <- compareCluster(ENTREZID~Cluster, data = data2,
                                fun="enrichKEGG", organism=org_code)
  p1 <- dotplot(formula_res, color ="qvalue", font.size = 7)
  keggenrich_name <-  paste0(dir_name, "/kegg_enrich.txt")
  write.table(as.data.frame(formula_res), file = keggenrich_name,
              row.names = F, col.names = T, sep = "\t", quote = F)

  formula_res_go <- compareCluster(ENTREZID~Cluster,
                                   data=data2, fun="enrichGO", OrgDb=org)
  g1 <- dotplot(formula_res_go, color ="qvalue", font.size = 7)

  goenrich_name <- paste0(dir_name, "/go_enrich.txt")
  write.table(as.data.frame(formula_res_go), file = goenrich_name,
              row.names = F, col.names = T, sep = "\t", quote = F)

  heatmap.file <- paste0(dir_name, "/heatmap.pdf")
  pdf(heatmap.file, width = 3, height = 4)
  print(ht)
  dev.off()
  kegg.file <- paste0(dir_name, "/enrichment_kegg.pdf")

  pdf(kegg.file, width = (km + 5/km) +1, height = 4)
  print(p1)
  dev.off()
  go.file <- paste0(dir_name, "/enrichment_go.pdf")
  pdf(go.file, width = (km + 12/km) +1, height = 4)
  print(g1)
  dev.off()
}



#' Gene symbol conversion from human to mouse
#'
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getLDS
#' @importFrom cowplot plot_grid
#' @param Gene_set Directory of Gene_set txt file
#' @export
#'
GeneSetConversion <- function(Gene_set) {
  dir_name <- Gene_set
  dir_name_1 <- paste(dir_name, "_convert", sep = "")
  dir.create(dir_name_1, showWarnings = F)
  gene_set_files <- list.files(path = Gene_set,
                               pattern = "*.txt")
  gene_set_files <- gsub(".txt", "", gene_set_files)
  gene_files_full <- list.files(path = Gene_set,
                                pattern = "*.txt", full.names = T)
  gene_files_full <- gsub(".txt", "", gene_files_full)
  gene_set_dir <- gsub(gene_set_files[1], "", gene_files_full[1])
  gene_files_full <- gsub(".txt", "", gene_files_full)

  for (name in gene_files_full) {
    data.file <- paste(name, ".txt", sep = "")
    print(data.file)
    genes <- read.table(data.file, skip = 2)
    genes <- genes$V1
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="asia.ensembl.org")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="asia.ensembl.org")
    genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                   values = genes ,mart = mouse,
                   attributesL = c("hgnc_symbol"),
                   martL = human, uniqueRows=T)
    gene.file <- gsub(gene_set_dir, "", data.file)
    gene.file2 <- paste(paste(dir_name_1, "/", sep = ""), gene.file, sep = "")
    write.table(genes, file = gene.file2, sep="\t", quote=F, row.names = F)
    genes <- read.table(gene.file2)
    set_name <- gsub(".txt", "", gene.file)
    data <- rbind(set_name, genes)
    write.table(data, file = gene.file2, sep="\t", quote=F, row.names = F, col.names = F)
  }
}

