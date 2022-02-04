#' Integrated heatmap
#'
#' @importFrom genefilter genescale
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @examples library(rnaseqviewer)
#'
#' data(DEG_list)
#'dir.create("DEG_list")
#' write.table(DEG_list[1], file = "DEG_list/dataset1.txt", sep = "\t", quote = FALSE)
#' write.table(DEG_list[2], file = "DEG_list/dataset2.txt", sep = "\t", quote = FALSE)
#' vennd("DEG_list")
#'
#' dir.create("count_list")
#' data(Row_count_data)
#' data(Row_count_data2)
#' write.table(Row_count_data, file = "count_list/data1.txt", sep = "\t", quote = FALSE)
#' write.table(Row_count_data2, file = "count_list/data2.txt", sep = "\t", quote = FALSE)
#'
#' int_heatmap(Count_matrix_dir = "count_list",
#'             Gene_set = "DEG_list/group_lists/dataset1:dataset2.txt")
#'
#' @references R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#' @param Count_matrix_dir Directory including count matrix txt files
#' @param Gene_set gene set txt file
#' @param pre_zscoring z-scoring before integration of data set
#' @export
#'
int_heatmap <- function(Count_matrix_dir, Gene_set, pre_zscoring = T){
matrix_files_full <- list.files(path = Count_matrix_dir,
                               pattern = "*.txt", full.names = T)
matrix_files <- list.files(path = Count_matrix_dir,
                          pattern = "*.txt")
matrix_files_name <- gsub(".txt", "", matrix_files)
matrix_dir <- gsub(matrix_files[1], "",matrix_files_full[1])
matrix_files_full <- gsub(".txt", "", matrix_files_full)
matrix_list <- list()
matrix_z_list <- list()
dir_name <- paste0(matrix_dir, "/integrated_Heatmap/")
dir.create(dir_name, showWarnings = F)
gene_set <- read.table(Gene_set, header = T, row.names = 1, sep = "\t")
for (name in matrix_files_name) {
  data.file <- paste0(matrix_dir, paste0(name, ".txt"))
  print(data.file)
  matrix <- read.table(data.file, header = T, row.names = 1, sep = "\t")
  matrix_2 <- matrix
  colnames(matrix_2) <- paste0(paste0(name, "-") ,colnames(matrix))
  if(pre_zscoring == T){
  matrix_z <- genescale(matrix_2, axis = 1, method = "Z")
  matrix_z <- na.omit(matrix_z)
  matrix_z <- merge(matrix, matrix_z, by = 0)[,-2:-(1 + length(colnames(matrix)))]
  matrix_z_list[name] <- list(matrix_z)
  }
  matrix_3 <- merge(matrix, matrix_2, by = 0)[,-2:-(1 + length(colnames(matrix)))]
  matrix_list[name] <- list(matrix_3)
}
base <- matrix_list[[1]]
int_matrix <- lapply(matrix_list[-1], function(i) base <<- merge(base, i, by = "Row.names"))
rownames(base) <- base$Row.names
base <- data.matrix(base[,-1])
pre_table.file <- paste0(paste0(dir_name, paste0(names(matrix_list), collapse = "-")),
                     "_ALL.txt")
write.table(base, file = pre_table.file, quote = F,sep = "\t", row.names = T, col.names = T)
base_2 <- merge(gene_set, base, by = 0)
if(length(colnames(gene_set)) != 0){
base_2 <- base_2[,-2:-(1 + length(colnames(gene_set)))]
}
rownames(base_2) <- base_2$Row.names
base_2 <- data.matrix(base_2[,-1])
if(pre_zscoring == T){
  base_z <- matrix_z_list[[1]]
  int_matrix <- lapply(matrix_z_list[-1], function(i) base_z <<- merge(base_z, i, by = "Row.names"))
  rownames(base_z) <- base_z$Row.names
  base_z <- data.matrix(base_z[,-1])
  all_z_table.file <- paste0(paste0(dir_name, paste0(names(matrix_list), collapse = "-")),
                       "_ALL_zscored.txt")
  write.table(base_z, file = all_z_table.file, quote = F,sep = "\t", row.names = T, col.names = T)
  base_z <- merge(gene_set, base_z, by = 0)
  if(length(colnames(gene_set)) != 0){
    base_z <- base_z[,-2:-(1 + length(colnames(gene_set)))]
  }
  rownames(base_z) <- base_z$Row.names
  base_z <- data.matrix(base_z[,-1])
}
if(pre_zscoring != T){
  base_2 <- genescale(base_2, axis = 1, method = "Z")
  base_z <- na.omit(base_2)
}
cond <- gsub(".+\\-", "", colnames(base_z))
cond <- gsub("\\_.+$", "", cond)
cond <- factor(cond, levels = unique(cond), ordered = TRUE)
if(length(rownames(base_z)) <= 50){
ht <- Heatmap(base_z, name = "z-score",
              clustering_method_columns = 'ward.D2',
              cluster_row_slices = T, show_row_names = T,
              top_annotation = HeatmapAnnotation(condition = cond))
}else{
  ht <- Heatmap(base_z, name = "z-score",
                clustering_method_columns = 'ward.D2',
                cluster_row_slices = T, show_row_names = F,
                top_annotation = HeatmapAnnotation(condition = cond))
}
Gene_set_name <- gsub(".+\\/","", Gene_set)
heatmap.file <- paste0(paste0(dir_name, paste0(names(matrix_list), collapse = "-")),
                       paste0(paste0("_", gsub(".txt", "", Gene_set_name)), ".pdf"))
pdf(heatmap.file,width = 7,height = 10)
print(ht)
dev.off()
pre_table.file2 <- paste0(paste0(dir_name, paste0(names(matrix_list), collapse = "-")),
                     paste0(paste0("_", gsub(".txt", "", Gene_set_name)), ".txt"))
write.table(base_2, file = pre_table.file2, quote = F,sep = "\t", row.names = T, col.names = T)
table.file <- paste0(paste0(dir_name, paste0(names(matrix_list), collapse = "-")),
                     paste0(paste0("_", gsub(".txt", "", Gene_set_name)), "_zscored.txt"))
write.table(base_z, file = table.file, quote = F,sep = "\t", row.names = T, col.names = T)

}
