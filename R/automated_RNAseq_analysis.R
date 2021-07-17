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
#' @param Count_matrix count matrix
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
    data <- data %>% tidyr::gather(key = sample,
                                   value = value, -Row.names)
    data$sample <- gsub("\\_.+$", "", data$sample)
    data$Row.names <- as.factor(data$Row.names)
    data$sample <- as.factor(data$sample)
    data$value <- as.numeric(data$value)
    stat.test <- data %>%
      group_by(Row.names) %>%
      tukey_hsd(value ~ sample) %>%
      add_significance("p.adj")
    stat.test <- stat.test %>% add_xy_position()
    stat.test
    image.file <- paste(name, ".pdf", sep = "")
    image.file <- gsub(group_dir, "", image.file)
    image.file <- paste(paste(dir_name_1, "/", sep = ""), image.file, sep = "")
    pdf(image.file, width = 15, height = 20)
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
    pdf(image.file2, width = 15, height = 20)
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
#' @param Count_matrix count matrix
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
#' @param Count_matrix Directory of count matrix
#' @param EBseq_Result result of EBseq
#' @param Species Species
#' @param Cond_1 Sumple number of condition_1
#' @param Cond_2 Sumple number of condition_2
#' @export
#'
pairwiseEBseq_viewer <- function(Count_matrix, EBseq_Result, Species, Cond_1, Cond_2){
  data <- read.table(EBseq_Result,header = T, row.names = 1)
  count <- read.table(Count_matrix,header=T, row.names = 1)
  data <- merge(data,count, by=0)

  switch (Species,
          "mouse" = org <- org.Mm.eg.db,
          "human" = org <- org.Hs.eg.db)
  switch (Species,
          "mouse" = org_code <- "mmu",
          "human" = org_code <- "hsa")

  my.symbols <- data$Row.names
  gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("ENTREZID", "SYMBOL", "UNIPROT"))
  colnames(gene_IDs) <- c("Row.names","ENTREZID","UNIPROT")
  data <- merge(data, gene_IDs, by="Row.names")
  data <- data %>% distinct(Row.names, .keep_all = T)
  log2FoldChange <- -1 * log2(data$PostFC)
  data <- cbind(data, log2FoldChange)
  baseMean <- (data$C1Mean + data$C2Mean)*(1/2)
  data <- cbind(data, baseMean)
  padj <- data$PPEE
  data <- cbind(data, padj)

  ##MA-plot
  m1 <- as.grob(ggmaplot(data, fdr = 0.05, fc = 2, size = 0.4,
                         palette = c("#B31B21", "#1465AC", "darkgray"),
                         genenames = as.vector(data$Row.names),
                         legend = "top", top = 20,
                         font.label = c("bold", 6),font.legend = "bold",
                         font.main = "bold",
                         ggtheme = ggplot2::theme_minimal(),
                         select.top.method = "fc"))
  data2 <- data[abs(data$log2FoldChange) > 1,]
  data2$group <- "upregulated"
  data2$group[data2$log2FoldChange < 0] <- "downregulated"
  data3 <- data2[abs(data2$PPEE) < 0.05,]

  ##heatmap
  data.z <- genescale(data3[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
  ht <- as.grob(Heatmap(data.z, name = "z-score",
                        clustering_method_columns = 'ward.D2',
                        show_row_names = F, show_row_dend = T))
  ma_name <- gsub("\\..+$", "", Count_matrix)
  ma_name <- paste(ma_name, "_ma-heatmap.pdf", sep = "")
  pdf(ma_name, width = 7, height = 3.5)
  print(plot_grid(m1, ht, rel_widths = c(2, 1)))
  dev.off()

  #dotplot
  universe <- AnnotationDbi::select(org,keys = rownames(count),
                                    keytype = "SYMBOL",
                                    columns = c("ENTREZID", "SYMBOL", "UNIPROT"))
  formula_res <- compareCluster(ENTREZID~group, data=data3,fun="enrichKEGG", organism=org_code, universe = universe)
  p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 9))
  keggenrich_name <- gsub("\\..+$", "", Count_matrix)
  keggenrich_name <- paste(keggenrich_name, "_kegg_enrich.csv", sep = "")
  write.table(as.data.frame(formula_res), file = keggenrich_name, row.names = F, col.names = T, sep = ",", quote = F)

  ##cnetplot
  upgene <- data3[data3$log2FoldChange > 1,]
  geneList_up <- upgene$log2FoldChange
  names(geneList_up) = as.character(upgene$ENTREZID)
  kk1 <- enrichKEGG(upgene$ENTREZID, organism =org_code,
                    pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  cnet1 <- setReadable(kk1, org, 'ENTREZID')
  p2 <- as.grob(cnetplot(cnet1, foldChange=geneList_up,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE))
  downgene <- data3[data3$log2FoldChange < 1,]
  geneList_down <- downgene$log2FoldChange
  names(geneList_down) = as.character(downgene$ENTREZID)
  kk2 <- enrichKEGG(downgene$ENTREZID, organism =org_code,
                    pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  cnet2 <- setReadable(kk2, org, 'ENTREZID')
  p3 <- as.grob(cnetplot(cnet2, foldChange=geneList_down,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE))

  ##GSEA plot
  data <- na.omit(data)
  geneList <- data$log2FoldChange
  names(geneList) = as.character(data$UNIPROT)
  geneList <- sort(geneList, decreasing = TRUE)
  barplot(sort(geneList, decreasing = T))
  kk2 <- gseKEGG(geneList = geneList, pvalueCutoff = 0.05,
                 organism = org_code, keyType = "uniprot",
                 exponent = 1, eps = 0, pAdjustMethod = "none",
                 minGSSize = 50, maxGSSize = 500, by = "fgsea",
                 use_internal_data = FALSE)
  p4 <- as.grob(gseaplot2(kk2, 1:6, pvalue_table = F))
  gsekegg_name <- gsub("\\..+$", "", Count_matrix)
  gsekegg_name <- paste(gsekegg_name, "_gsekegg.csv", sep = "")
  write.table(as.data.frame(kk2), file = gsekegg_name, row.names = F, col.names = T, sep = ",", quote = F)

  kegg_name <- gsub("\\..+$", "", Count_matrix)
  kegg_name <- paste(gsekegg_name, "_kegg.pdf", sep = "")
  pdf(kegg_name, width = 14, height = 14.0)
  print(plot_grid(p1, p4, p2, p3, nrow = 2))
  dev.off()



  ##GOenrichment
  ##enrichment
  formula_res_go <- compareCluster(ENTREZID~group,
                                   data=data3,fun="enrichGO", OrgDb=org)
  g1 <- as.grob(dotplot(formula_res_go, color ="qvalue", font.size = 9))
  goenrich_name <- gsub("\\..+$", "", Count_matrix)
  goenrich_name <- paste(goenrich_name, "_go_enrich.csv", sep = "")
  write.table(as.data.frame(formula_res_go), file = goenrich_name,
              row.names = F, col.names = T, sep = ",", quote = F)
  ##cnetplot
  go1 <- enrichGO(upgene$ENTREZID, OrgDb = org,
                  pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  cnet_go1 <- setReadable(go1, org, 'ENTREZID')
  g2 <- as.grob(cnetplot(cnet_go1, foldChange=geneList_up,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE))
  go2 <- enrichGO(downgene$ENTREZID, OrgDb = org,
                  pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
  cnet_go2 <- setReadable(go2, org, 'ENTREZID')
  g3 <- as.grob(cnetplot(cnet_go2, foldChange=geneList_down,
                         cex_label_gene = 0.5, cex_label_category = 0.75,
                         cex_category = 0.5, colorEdge = TRUE))

  ##GSEA plot
  data <- na.omit(data)
  geneList <- data$log2FoldChange
  names(geneList) = as.character(data$ENTREZID)
  geneList <- sort(geneList, decreasing = TRUE)
  barplot(sort(geneList, decreasing = T))
  go3 <- gseGO(geneList = geneList, pvalueCutoff = 0.05,
               OrgDb = org, exponent = 1, eps = 0,
               pAdjustMethod = "none", minGSSize = 50,
               maxGSSize = 500, by = "fgsea")
  g4 <- as.grob(gseaplot2(kk2, 1:6, pvalue_table = F))
  gsego_name <- gsub("\\..+$", "", Count_matrix)
  gsego_name <- paste(gsego_name, "_gseGO.csv", sep = "")
  write.table(as.data.frame(kk2), file = gsego_name, row.names = F, col.names = T, sep = ",", quote = F)
  go_name <- gsub("\\..+$", "", Count_matrix)
  go_name <- paste(go_name, "_GO.pdf", sep = "")
  pdf(go_name, width = 14, height = 14.0)
  print(plot_grid(g1, g4, g2, g3, nrow = 2))
  dev.off()

  #FC上位50下位50をboxplot
  data4 <- data3[sort(data3$log2FoldChange, decreasing = T, index=T)$ix,]
  up50 <- data4[1:50,8:(7 + Cond_1 + Cond_2)]
  up50$Row.names <- data4[1:50,]$Row.names
  up50 <- up50 %>% gather(key=sample, value=value,-Row.names)
  up50$sample <- gsub("\\_.+$", "", up50$sample)
  up50$Row.names <- as.factor(up50$Row.names)
  up50$sample <- as.factor(up50$sample)
  up50$value <- as.numeric(up50$value)
  data4 <- data3[sort(data3$log2FoldChange, decreasing = F, index=T)$ix,]
  down50 <- data4[1:50,8:(7 + Cond_1 + Cond_2)]
  down50$Row.names <- data4[1:50,]$Row.names
  down50 <- down50 %>% gather(key=sample, value=value,-Row.names)
  down50$sample <- gsub("\\_.+$", "", down50$sample)
  down50$Row.names <- as.factor(down50$Row.names)
  down50$sample <- as.factor(down50$sample)
  down50$value <- as.numeric(down50$value)

  data5 <- data4[,8:(7 + Cond_1 + Cond_2)]
  rownames(data5) <- data4$Row.names
  deg_name <- gsub("\\..+$", "", Count_matrix)
  deg_name <- paste(deg_name, "_DEG_count.csv", sep = "")
  write.table(data5, file = deg_name, row.names = T, col.names = T, sep = ",", quote = F)

  degtop_name <- gsub("\\..+$", "", Count_matrix)
  degtop_name <- paste(degtop_name, "_top_DEG.pdf", sep = "")
  pdf(degtop_name,height = 10, width = 10)
  plot(ggpubr::ggboxplot(up50, x = "sample", y = "value",
                         fill = "sample", facet.by = "Row.names",
                         scales = "free", add = "jitter",
                         xlab = "gene", ylab = "TPM")+
         ggplot2::theme(axis.text.x= ggplot2::element_text(size = 5),
                        axis.text.y= ggplot2::element_text(size = 10)) +
         ggplot2::scale_y_continuous(limits = c(0, NA)))
  plot(ggpubr::ggboxplot(down50, x = "sample", y = "value",
                         fill = "sample", facet.by = "Row.names",
                         scales = "free", add = "jitter",
                         xlab = "gene", ylab = "TPM")+
         ggplot2::theme(axis.text.x= ggplot2::element_text(size = 5),
                        axis.text.y= ggplot2::element_text(size = 10)) +
         ggplot2::scale_y_continuous(limits = c(0, NA)))
  dev.off()
}


