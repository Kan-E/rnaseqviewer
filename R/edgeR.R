#' edgeR
#'
#' @importFrom edgeR DGEList
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR estimateCommonDisp
#' @importFrom edgeR estimateTagwiseDisp
#' @importFrom edgeR exactTest
#' @importFrom edgeR topTags
#' @importFrom edgeR filterByExpr
#' @importFrom qvalue qvalue
#' @importFrom IHW ihw
#' @importFrom IHW as.data.frame
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @examples library(rnaseqviewer)
#'
#' data(Row_count_data)
#' write.table(Row_count_data, file = "Row_count_data.txt", sep = "\t", quote = FALSE)
#' edger("Row_count_data.txt")
#'
#' @references Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#' @references Nikolaos Ignatiadis, Bernd Klaus, Judith Zaugg and Wolfgang Huber (2016): Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. Nature Methods 13:577, doi: 10.1038/nmeth.3885
#' @references John D. Storey, Andrew J. Bass, Alan Dabney and David Robinson (2021). qvalue: Q-value estimation for false discovery rate control. R package version 2.26.0. http://github.com/jdstorey/qvalue
#' @param Row_count_matrix Row Count matrix txt file (Not normalized count matrix)
#' @param method BH, Qvalue, or IHW
#' @export
#'
edger <- function(Row_count_matrix, method = "BH"){
  count<-read.table(Row_count_matrix, header = T, row.names = 1, sep = "\t")
  count <- as.matrix(count)
  if(length(grep("/", Row_count_matrix) == 1)){
    dir_name <- paste0(gsub("/[^/]+$", "", Row_count_matrix), "/")
    file_name <- gsub(gsub("/[^/]+$", "", Row_count_matrix), "", Row_count_matrix)
    file_name <- gsub("/", "",file_name)
  }else{
    dir_name <- ""
    file_name <- Row_count_matrix
  }
  file_name <- gsub("\\..+$", "", file_name)
  collist <- gsub("\\_.+$", "", colnames(count))
  group <- factor(collist)
  name <- paste0(paste0(unique(collist)[1], "-vs-"), paste0(unique(collist)[2],  paste0("_edgeR-",method)))
  print(name)
  d <- DGEList(counts = count, group = group)
  keep <- filterByExpr(d)
  d = d[keep, , keep.lib.sizes=FALSE]
  d <- calcNormFactors(d)
  d <- estimateCommonDisp(d)
  d <- estimateTagwiseDisp(d)
  norm_counts.table <- t(t(d$pseudo.counts)*(d$samples$norm.factors))
  count_file <- paste0(paste0(dir_name, paste0("Normalized_count_matrix_from_"), file_name),
                        paste0(paste0("_", name), ".txt"))
  write.table(norm_counts.table, file=count_file, sep="\t", quote=F, row.names = T)
  result <- exactTest(d, pair = c(unique(group)[2],unique(group)[1]))
  table <- as.data.frame(topTags(result, n = nrow(count)))
  qvalue <- qvalue::qvalue(table$PValue)
  table$padj <- qvalue$qvalues
  ihw_res <- ihw(PValue ~ 2^logCPM,  data=table, alpha = 0.1)
  ihw_res_df <- IHW::as.data.frame(ihw_res)
  table$ihw_padj <- ihw_res_df$adj_pvalue
  if(method == "BH"){label <- c("log2FoldChange", "log2CPM", "PValue","padj", "Qvalue", "IHW_FDR")}
  if(method == "Q"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "padj", "IHW_FDR")}
  if(method == "IHW"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "Qvalue", "padj")}
  colnames(table) <- label
  result_file <- paste0(paste0(dir_name, paste0("result_of_"), file_name),
                        paste0(paste0("_", name), ".txt"))
  write.table(table, file = result_file, row.names = T, col.names = T, sep = "\t", quote = F)
}
