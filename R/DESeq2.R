#' deseq2
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom DESeq2 counts
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @examples library(rnaseqviewer)
#'
#' data(Row_count_data)
#' write.table(Row_count_data, file = "Row_count_data.txt", sep = "\t")
#' deseq2("Row_count_data.txt")
#'
#' @docType data
#' @name Row_count_data
#' @references Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#' @param Row_count_matrix Row Count matrix txt file (Not normalized count matrix)
#' @export
#'
deseq2 <- function(Row_count_matrix){
count<-read.table(Row_count_matrix, header = T, row.names = 1)
collist <- gsub("\\_.+$", "", colnames(count))
group <- data.frame(con = factor(collist))
name <- paste0(paste0(unique(collist)[1], "-vs-"), paste0(unique(collist)[2], "_DEseq2"))
print(name)
dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
dds <- DESeq(dds)
res <- results(dds)
result_file <- paste0("result_of_", paste0(name, ".txt"))
write.table(res, file = result_file, row.names = T, col.names = T, sep = "\t", quote = F)
normalized_counts <- counts(dds, normalized=TRUE)
count_file <- paste0("Normalized_count_matrix_from_", paste0(name, ".txt"))
write.table(normalized_counts, file = count_file, row.names = T, col.names = T, sep = "\t", quote = F)
}
