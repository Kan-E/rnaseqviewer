#' deseq2
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom DESeq2 counts
#' @importFrom IHW ihw
#' @importFrom IHW as.data.frame
#' @importFrom qvalue qvalue
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @examples library(rnaseqviewer)
#'
#' data(Row_count_data)
#' write.table(Row_count_data, file = "Row_count_data.txt", sep = "\t", quote = FALSE)
#' deseq2("Row_count_data.txt")
#'
#' @references Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' @references Nikolaos Ignatiadis, Bernd Klaus, Judith Zaugg and Wolfgang Huber (2016): Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. Nature Methods 13:577, doi: 10.1038/nmeth.3885
#' @references John D. Storey, Andrew J. Bass, Alan Dabney and David Robinson (2021). qvalue: Q-value estimation for false discovery rate control. R package version 2.26.0. http://github.com/jdstorey/qvalue
#' @param Row_count_matrix Row Count matrix txt file (Not normalized count matrix)
#' @param method BH, Qvalue, or IHW
#' @export
#'
deseq2 <- function(Row_count_matrix, method = "BH"){
count<-read.table(Row_count_matrix, header = T, row.names = 1, sep = "\t")

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
group <- data.frame(con = factor(collist))
name <- paste0(paste0(unique(collist)[1], "-vs-"), paste0(unique(collist)[2], paste0("_DEseq2-",method)))
print(name)
dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
dds$con <- factor(dds$con, levels = unique(collist))
dds <- DESeq(dds)
contrast <- c("con", unique(collist))
res <- results(dds,  contrast = contrast)
if(method == "IHW") {
  ihw_res <- ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha = 0.1)
  res$padj <- IHW::as.data.frame(ihw_res)$adj_pvalue
}
if(method == "Qvalue") {
  res <- results(dds,  contrast = contrast)
  qvalue <- qvalue::qvalue(res$pvalue)
  res$padj <- qvalue$qvalues
}
result_file <- paste0(paste0(dir_name, paste0("result_of_"), file_name),
                      paste0(paste0("_", name), ".txt"))
write.table(res, file = result_file, row.names = T, col.names = T, sep = "\t", quote = F)
normalized_counts <- counts(dds, normalized=TRUE)
count_file <- paste0(paste0(dir_name, paste0("Normalized_count_matrix_from_"), file_name),
                     paste0(paste0("_", name), ".txt"))
write.table(normalized_counts, file = count_file, row.names = T, col.names = T, sep = "\t", quote = F)
}
