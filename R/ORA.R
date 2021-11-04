#' Over-representation analysis (ORA) by using GO and KEGG gene sets
#'
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
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
#' @importFrom cowplot plot_grid
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
#' ORA("DEG_list/group_lists", Species = "human")
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
#' @param gene_list_dir Directory including gene list txt files
#' @param Species Species
#' @param color qvalue or p.adjust
#' @export
#'
ORA <- function(gene_list_dir, Species, color = "qvalue") {
  data_files_full <- list.files(path = gene_list_dir,
                                pattern = "*.txt", full.names = T)
  data_files <- list.files(path = gene_list_dir,
                           pattern = "*.txt")
  data_dir <- gsub(data_files[1], "",data_files_full[1])
  currentD <- getwd()
  setwd(data_dir)
  dir.create("result", showWarnings = F)
  data_files <- gsub(".txt", "", data_files)
  df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
  for (name in data_files) {
    data.file <- paste(name, ".txt", sep = "")
    print(data.file)
    data <- read.table(data.file, header = T,
                       row.names = 1)
    df2 <- data.frame(GeneID = rownames(data), Group = name)
    df <- rbind(df, df2)
  }
  setwd(currentD)
  if (!is.null(Species)){
    switch (Species,
            "mouse" = org <- org.Mm.eg.db,
            "human" = org <- org.Hs.eg.db)
    switch (Species,
            "mouse" = org_code <- "mmu",
            "human" = org_code <- "hsa")
  }
    my.symbols <- df$GeneID
    gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                      keytype = "SYMBOL",
                                      columns = c("ENTREZID", "SYMBOL"))
    colnames(gene_IDs) <- c("GeneID","ENTREZID")
    data2 <- merge(df, gene_IDs, by="GeneID")
    formula_res <- try(compareCluster(ENTREZID~Group, data = data2,
                                  fun="enrichKEGG", organism=org_code), silent = T)
    if (class(formula_res) == "try-error") {
      formula_res <- NA
      p1 <- NULL
    } else {
      if ((length(as.data.frame(formula_res)) == 0) ||
          is.na(unique(as.data.frame(formula_res)$qvalue))) {
        p1 <- NULL
      } else{
    p1 <- as.grob(dotplot(formula_res, font.size = 7, color = color))
    keggenrich_name <-  paste0(data_dir, "result/kegg_enrich.txt")
    write.table(as.data.frame(formula_res), file = keggenrich_name,
                row.names = F, col.names = T, sep = "\t", quote = F)
      }
    }
    cnetkegg_list <- list()
    for (name in unique(data2$Group)) {
      kk1 <- enrichKEGG(data2$ENTREZID[data2$Group == name], organism =org_code)
      if(is.null(kk1)){
        cnet1 <- NULL
      } else cnet1 <- setReadable(kk1, org, 'ENTREZID')
      if (length(cnet1$ID) == 0) {
        p2 <- NULL
      } else{
        p2 <- as.grob(cnetplot(cnet1, cex_label_gene = 0.5, cex_label_category = 0.75,showCategory = 5,
                               cex_category = 0.5, colorEdge = TRUE)+ guides(edge_color = "none"))
        cnetkegg_list[[name]] = p2
        keggcnet_name <-  paste0(data_dir, "result/kcnet_")
        keggcnet_name2 <-  paste0(keggcnet_name, paste0(name, ".pdf"))
        pdf(keggcnet_name2, width = 5, height = 5)
        print(plot_grid(p2))
        dev.off()
        kcnet_name <-  paste0(keggcnet_name, paste0(name, ".txt"))
        write.table(as.data.frame(cnet1), file = kcnet_name,
                    row.names = F, col.names = T, sep = "\t", quote = F)
      }
    }

    formula_res_go <- try(compareCluster(ENTREZID~Group,
                                     data=data2, fun="enrichGO", OrgDb=org),silent = T)
    if (class(formula_res_go) == "try-error") {
      formula_res_go <- NA
      g1 <- NULL
    } else{
      if ((length(as.data.frame(formula_res_go)) == 0) ||
          is.na(unique(as.data.frame(formula_res_go)$qvalue))) {
        g1 <- NULL
      } else{
    g1 <- as.grob(dotplot(formula_res_go, font.size = 7, color = color))
    goenrich_name <- paste0(data_dir, "/result/go_enrich.txt")
    write.table(as.data.frame(formula_res_go), file = goenrich_name,
                row.names = F, col.names = T, sep = "\t", quote = F)
      }
    }
    cnetgo_list <- list()
    for (name in unique(data2$Group)) {
      go1 <- enrichGO(data2$ENTREZID[data2$Group == name], OrgDb = org)
      if(is.null(go1)){
        cnet_go1 <- NULL
      } else cnet_go1 <- setReadable(go1, org, 'ENTREZID')
      if (length(cnet_go1$ID) == 0) {
        g2 <- NULL
      } else{
        g2 <- as.grob(cnetplot(cnet_go1, cex_label_gene = 0.5, cex_label_category = 0.75,
                               cex_category = 0.5, colorEdge = TRUE)+ guides(edge_color = "none"))
        cnetgo_list[[name]] = g2
        gocnet_name <-  paste0(data_dir, "result/gcnet_")
        gocnet_name2 <-  paste0(gocnet_name, paste0(name, ".pdf"))
        pdf(gocnet_name2, width = 5, height = 5)
        print(plot_grid(g2))
        dev.off()
        gcnet_name <-  paste0(gocnet_name, paste0(name, ".txt"))
        write.table(as.data.frame(cnet_go1), file = gcnet_name,
                    row.names = F, col.names = T, sep = "\t", quote = F)
      }
    }

    kegg.file <- paste0(data_dir, "/result/enrichment_kegg.pdf")
    n <- length(data_files)
    pdf(kegg.file, width = (n + 5/n) +1, height = 5)
    print(plot_grid(p1, ncol = 1))
    dev.off()
    go.file <- paste0(data_dir, "/result/enrichment_go.pdf")
    pdf(go.file, width = (n + 12/n) +1, height = 5)
    print(plot_grid(g1, ncol = 1))
    dev.off()
}
