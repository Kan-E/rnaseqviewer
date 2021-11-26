#' Visualization of multiple comparison DEG analysis
#'
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom rstatix group_by
#' @importFrom ggpubr ggboxplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 unit
#' @importFrom ggplot2 guides
#' @importFrom tidyr gather
#' @importFrom dplyr %>%
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom AnnotationDbi select
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ggplotify as.grob
#' @importFrom gridExtra grid.arrange
#' @importFrom gridExtra arrangeGrob
#' @importFrom clusterProfiler compareCluster
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler enrichGO
#' @importFrom enrichplot dotplot
#' @importFrom enrichplot cnetplot
#' @importFrom DOSE setReadable
#' @importFrom graphics barplot
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @import ggnewscale
#' @importFrom cowplot plot_grid
#' @examples library(rnaseqviewer)
#'
#' #' #three conditions DEG analysis
#'
#' data("Row_count_3conditions")
#' write.table(Row_count_3conditions, file = "Row_count_3conditions.txt", sep = "\t", quote = FALSE)
#' ebseq("Row_count_3conditions.txt")
#'
#' multiDEG_overview(Normalized_count_matrix =
#'                   "Normalized_count_matrix_from_Cond1-vs-Cond2-vs-Cond3_EBseq.txt",
#'                   EBseq_result = "result_of_Cond1-vs-Cond2-vs-Cond3_EBseq.txt",
#'                   EBseq_condmeans ="result_of_Cond1-vs-Cond2-vs-Cond3_EBseq.condmeans",
#'                   Species = "human", fdr = 0.05, fc = 1.5, basemean = 0.5)
#'
#' @references T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
#' @references Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609
#' @references Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0.
#' @references Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.
#' @references Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.
#' @references R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#' @references Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr
#' @references Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
#' @param Normalized_count_matrix Count matrix txt file (e.g. TPM count matrix.txt)
#' @param EBseq_result result txt file of EBseq analysis
#' @param EBseq_condmeans Condmeands txt file from EBseq analysis
#' @param Species Species
#' @param fdr Accepted false discovery rate for considering genes as differentially expressed
#' @param fc the fold change threshold. Only genes with a fold change >= fc and padj <= fdr are considered as significantly differentially expressed.
#' @param basemean basemean threshold.
#' @export
#'
multiDEG_overview <- function(Normalized_count_matrix, EBseq_result, EBseq_condmeans,
                           Species = NULL, fdr = 0.05, fc = 2, basemean = 0){
  data <-read.table(Normalized_count_matrix,header = T, sep = "\t")
  collist <- gsub("\\_.+$", "", colnames(data))
  vec <- c()
  for (i in 1:length(unique(collist))) {
    num <- length(collist[collist == unique(collist)[i]])
    vec <- c(vec, num)
  }
  Cond_1 <- vec[1]
  Cond_2 <- vec[2]
  Cond_3 <- vec[3]
  collist <- unique(collist[-1])
  result_Condm <- read.table(EBseq_condmeans, header = T, row.names = 1, sep = "\t")
  result_FDR <- read.table(EBseq_result,header = T, row.names = 1, sep = "\t")

  dir_name <- gsub(".txt", "", Normalized_count_matrix)
  dir_name <- paste0(dir_name, paste0("_fc", fc))
  dir_name <- paste0(dir_name, paste0("_fdr", fdr))
  dir_name <- paste0(dir_name, paste0("_basemean", basemean))
  dir.create(dir_name, showWarnings = F)

  if (!is.null(Species)){
  switch (Species,
          "mouse" = org <- org.Mm.eg.db,
          "human" = org <- org.Hs.eg.db)
  switch (Species,
          "mouse" = org_code <- "mmu",
          "human" = org_code <- "hsa")
  }
  plot_list = list()
  dotplot_list = list()
  htplot_list = list()
  boxplot_list = list()
  data4_sum = data.frame()
  for (i in 1:3) {
    if(i == 1) {specific = collist[1]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
    result_Condm$FC_x <- log2((result_Condm$C1 + 0.01)/(result_Condm$C2 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C1 + 0.01)/(result_Condm$C3 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    }
    if(i == 2) {specific = collist[2]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
    result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
    Pattern1 <- "Pattern3"
    Pattern2 <- "Pattern5"
    }
    if(i == 3) {specific = collist[3]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
    result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
    Pattern1 <- "Pattern2"
    Pattern2 <- "Pattern5"
    }
    dir_name2 <- paste0(dir_name, paste0("/", specific))
    dir.create(dir_name2, showWarnings = F)
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > basemean)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= fdr & result$FC_x < log2(1/fc) & result$FC_y < log2(1/fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= fdr & result$FC_x > log2(fc) & result$FC_y > log2(fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
    data3 <- data.frame(Row.names = result$Row.names, FC_x = result$FC_x,
                        FC_y = result$FC_y, padj = result$FDR, sig = sig, FC_xy = result$FC_x * result$FC_y)
    if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
      new.levels <- c( paste0(paste0(specific,"_up: "), sum(sig == 1)), paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("red","blue", "darkgray")}
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      new.levels <- c(paste0(paste0(specific,"_up: "), sum(sig == 1)), "NS" )
      col = c("red", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      new.levels <- c(paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("blue", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) == 0)){
      new.levels <- c("NS")
      col = "darkgray"}

    data3$sig <- factor(data3$sig, labels = new.levels)
    complete_data <- stats::na.omit(data3)
    labs_data <- subset(complete_data, padj <= fdr & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(fc))
    labs_data<-  labs_data[sort(labs_data$FC_xy, decreasing = T, index=T)$ix,]
    labs_data <- dplyr::filter(labs_data, sig != "NS")
    labs_data <- utils::head(labs_data, 20)
    font.label <- data.frame(size=5, color="black", face = "plain")
    set.seed(42)
    FC_x <- FC_y <- sig <- Row.names <- padj <- NULL
    p <- ggplot(data3, aes(x = FC_x, y = FC_y)) + geom_point(aes(color = sig),size = 0.1)
    p <- p  + geom_hline(yintercept = c(-log2(fc), log2(fc)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(fc), log2(fc)),linetype = c(2, 2), color = c("black", "black"))
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = Row.names),
                                      box.padding = unit(0.35, "lines"), point.padding = unit(0.3,
                                                                                              "lines"), force = 1, fontface = font.label$face,
                                      size = font.label$size/3, color = font.label$color)
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 7),
            axis.text.y= ggplot2::element_text(size = 7),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 8)) +
      xlab(FC_xlab) + ylab(FC_ylab)
    plot_list[[i]] = p

    if(length(new.levels) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      data4_sum = rbind(data4_sum, data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if (!is.null(Species)){
      my.symbols <- data2$Row.names
      SYMBOL <- NULL
      gene_IDs<-AnnotationDbi::select(org, keys = my.symbols,keytype = "SYMBOL",columns = c("ENTREZID", "SYMBOL"))
      colnames(gene_IDs) <- c("Row.names","ENTREZID")
      data4 <- merge(data4, gene_IDs, by="Row.names")
      data4 <- data4 %>% dplyr::distinct(Row.names, .keep_all = T)
      data4_sum = rbind(data4_sum, data4)
      universe <- AnnotationDbi::select(org, keys = rownames(data),keytype = "SYMBOL",columns = c("ENTREZID", "SYMBOL"))
      universe <- universe %>% distinct(SYMBOL, .keep_all = T)

      formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichKEGG", organism=org_code, universe = universe), silent = T)
      if (class(formula_res) == "try-error") {
        formula_res <- NA
        d <- NULL
      } else {
        if ((length(as.data.frame(formula_res)) == 0) ||
            is.na(unique(as.data.frame(formula_res)$qvalue))) {
          d <- NULL
        } else{
          d <- as.grob(dotplot(formula_res, showCategory=5, color ="qvalue" ,font.size=7))
          keggenrich_name <- paste0(paste0(dir_name2, "/"),
                                    paste0(specific,"_sig_keggenrich.txt"))
          write.table(as.data.frame(formula_res), file = keggenrich_name,
                      row.names = F, col.names = T, sep = "\t", quote = F)
        }
      }
      formula_res_go <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichGO", OrgDb=org, universe = universe), silent = T)
      if (class(formula_res_go) == "try-error") {
        formula_res_go <- NA
        g1 <- NULL
      } else {
        if ((length(as.data.frame(formula_res_go)) == 0) ||
                  is.na(unique(as.data.frame(formula_res_go)$qvalue))) {
        g1 <- NULL
      } else{
        g1 <- as.grob(dotplot(formula_res_go, color ="qvalue", font.size = 7))
        goenrich_name <- paste0(paste0(dir_name2, "/"),
                                paste0(specific,"_sig_goenrich.txt"))
        write.table(as.data.frame(formula_res_go), file = goenrich_name,
                    row.names = F, col.names = T, sep = "\t", quote = F)
      }
      }

      cnetkegg_list <- list()
      cnetgo_list <- list()
    }
      for (name in new.levels) {
        if (name != "NS"){
          if (!is.null(Species)){
          kk1 <- enrichKEGG(data4$ENTREZID[data4$sig == name], organism =org_code)
          if (is.null(kk1)) {
            cnet1 <- NULL
          } else cnet1 <- setReadable(kk1, org, 'ENTREZID')
          if ((length(cnet1$ID) == 0) || is.na(unique(cnet1$qvalue))) {
            c <- NULL
          } else{
            c <- cnetplot(cnet1, cex_label_gene = 0.5, cex_label_category = 0.75,
                          cex_category = 0.5, colorEdge = TRUE)
            c <- as.grob(c + guides(edge_color = "none"))
            cnetkegg_list[[name]] = c
          }
          go1 <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb = org)
          if (is.null(go1)) {
            cnet_go1 <- NULL
          } else cnet_go1 <- setReadable(go1, org, 'ENTREZID')
          if ((length(cnet_go1$ID) == 0) || is.na(unique(cnet_go1$qvalue))) {
            g <- NA
          } else{
            g <- cnetplot(cnet_go1, cex_label_gene = 0.5, cex_label_category = 0.75,
                          cex_category = 0.5, colorEdge = TRUE)
            g <- as.grob(g + guides(edge_color = "none"))
            cnetgo_list[[name]] = g
          }
          }

          ##boxplot
          boxdata <- dplyr::filter(data4, sig == name)
          boxdata <- data.frame(Row.names = boxdata$Row.names, FC_xy = boxdata$FC_xy)
          boxdata2 <- merge(boxdata, data, by = "Row.names")
          rownames(boxdata2) <- boxdata2$Row.names
          boxdata3 <- boxdata2[sort(boxdata2$FC_xy, decreasing = T, index=T)$ix,]
          allup <- boxdata3[,-2]
          up50 <- boxdata3[1:50, -2]
          collist <- gsub("\\_.+$", "", colnames(up50))
          collist <- unique(collist[-1])
          rowlist <- gsub("\\_.+$", "", up50[,1])
          rowlist <- unique(rowlist)
          value <- NULL
          up50 <- up50 %>% gather(key=sample, value=value,-Row.names)
          up50$sample <- gsub("\\_.+$", "", up50$sample)
          up50$Row.names <- as.factor(up50$Row.names)
          up50$value <- as.numeric(up50$value)
          up50$sample <- factor(up50$sample,levels=collist,ordered=TRUE)
          up <- ggpubr::ggboxplot(up50,x = "sample", y = "value",fill = "sample",
                                  facet.by = "Row.names",scales = "free",
                                  add = "jitter", add.params = list(size=0.5),
                                  xlab = FALSE, legend = "none", ylim = c(0, NA))+
            ggplot2::theme(axis.text.x= ggplot2::element_text(size = 5),
                           axis.text.y= ggplot2::element_text(size = 7),
                           panel.background = element_rect(size = 0.5),
                           title = ggplot2::element_text(size = 7),
                           text = ggplot2::element_text(size = 10)) +
            scale_fill_manual(values=c("gray", "#4dc4ff", "#ff8082"))
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
          if (length(rowlist) <= 2) pdf_size <- 3
          if (length(rowlist) > 100) pdf_size <- 16.5
          top50_boxplot_name <- paste0(paste0(dir_name2, "/"),
                                       paste0(name,"_top50_boxplot.pdf"))
          pdf(top50_boxplot_name, height = pdf_size, width = pdf_size)
          suppressWarnings(print(up))
          dev.off()
          table_name <- paste0(paste0(dir_name2, "/"),
                               paste0(name, " genes_count.txt"))
          write.table(allup, file = table_name, row.names = F, col.names = T, sep = "\t", quote = F)
        }
      }
      if (!is.null(Species)){
      if (length(cnetkegg_list) == 2){
        cnetkegg1 <- cnetkegg_list[[1]]
        cnetkegg2 <- cnetkegg_list[[2]]}
      if (length(cnetkegg_list) == 1){
        cnetkegg1 <- cnetkegg_list[[1]]
        cnetkegg2 <- NULL}
      if (length(cnetkegg_list) == 0){
        cnetkegg1 <- NULL
        cnetkegg2 <- NULL}
      if (length(cnetgo_list) == 2){
        cnetgo1 <- cnetgo_list[[1]]
        cnetgo2 <- cnetgo_list[[2]]}
      if (length(cnetgo_list) == 1){
        cnetgo1 <- cnetgo_list[[1]]
        cnetgo2 <- NULL}
      if (length(cnetgo_list) == 0){
        cnetgo1 <- NULL
        cnetgo2 <- NULL}
      keggenrichplot_name <- paste0(paste0(dir_name2, "/"),
                                    paste0(specific,"_sig_enrichplot.pdf"))
      pdf(keggenrichplot_name, height = 8, width = 15)
      print(plot_grid(d, cnetkegg1, cnetkegg2,
                      g1, cnetgo1, cnetgo2, ncol =3, nrow = 2))
      dev.off()
      }

      data5 <- merge(data4, data, by ="Row.names")
      rownames(data5) <- data5$Row.names
      if (!is.null(Species)){
        data5 <- data5[,-1:-7]
      } else {data5 <- data5[,-1:-6]}
      data.z <- genescale(data5, axis=1, method="Z")
      ht <- as.grob(Heatmap(data.z, name = "z-score", column_order = colnames(data.z),
                            clustering_method_columns = 'ward.D2',
                            show_row_names = F, show_row_dend = T))
      htplot_list[[i]] = ht
      table_name2 <- paste0(paste0(dir_name2, "/"),
                           paste0(specific, " genes_sig_count.txt"))
      write.table(data5, file = table_name2, row.names = T, col.names = T, sep = "\t", quote = F)
    }
  }
  scatterplot_heatmap_name <- paste0(paste0(dir_name, "/"), "scatterplot_heatmap.pdf")
  pdf(scatterplot_heatmap_name, height = 10.5, width = 6)
  suppressWarnings(print(plot_grid(plot_list[[1]], htplot_list[[1]],
                  plot_list[[2]], htplot_list[[2]],
                  plot_list[[3]], htplot_list[[3]],
                  nrow = 3, ncol = 2)))
  dev.off()
}
