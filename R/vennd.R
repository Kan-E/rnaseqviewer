#' Venn Diagram
#'
#' @importFrom venn venn
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @param gene_list_dir Directory including gene list txt files
#' @references Adrian Dusa (2021). venn: Draw Venn Diagrams. R package version 1.10. https://CRAN.R-project.org/package=venn
#' @export
#'
vennd <- function(gene_list_dir) {
  data_files_full <- list.files(path = gene_list_dir,
                                pattern = "*.txt", full.names = T)
  data_files <- list.files(path = gene_list_dir,
                           pattern = "*.txt")
  data_dir <- gsub(data_files[1], "",data_files_full[1])
  currentD <- getwd()
  setwd(data_dir)
  dir.create("group_lists", showWarnings = F)
  data_files <- gsub(".txt", "", data_files)
  gene_list = list()
  for (name in data_files) {
    data.file <- paste(name, ".txt", sep = "")
    print(data.file)
    data <- read.table(data.file, header = T,
                       row.names = 1)
    gene_list[name] <- list(rownames(data))
  }
  setwd(currentD)
  venn.file <- paste0(data_dir, "venn.pdf")
  pdf(venn.file, height = 3, width = 3)
  venn(gene_list, ilab=TRUE, zcolor = "style", ilcs = 0.8, sncs = 0.6 )
  dev.off()

  df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
  for (name in names(attr(venn(gene_list),"intersections"))){
    data3 <- data.frame(Row.names = attr(venn(gene_list),"intersections")[name])
    table.file <- paste0(paste0(data_dir, "group_lists/"), paste0(name, ".txt"))
    write.table(data3, file = table.file, row.names = F, col.names = T, sep = "\t", quote = F)
  }
 }
