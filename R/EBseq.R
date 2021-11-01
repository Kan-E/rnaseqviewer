#' ebseq
#'
#' @importFrom EBSeq MedianNorm
#' @importFrom EBSeq GetNormalizedMat
#' @importFrom EBSeq EBTest
#' @importFrom EBSeq PostFC
#' @importFrom EBSeq GetPatterns
#' @importFrom EBSeq EBMultiTest
#' @importFrom EBSeq GetMultiPP
#' @importFrom EBSeq GetMultiFC
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @examples library(rnaseqviewer)
#'
#' #pairwise DEG analysis
#'
#' data(Row_count_data)
#' write.table(Row_count_data, file = "Row_count_data.txt", sep = "\t")
#' ebseq("Row_count_data.txt")
#'
#' #three conditions DEG analysis
#'
#' data("Row_count_3conditions")
#' write.table(Row_count_3conditions, file = "Row_count_3conditions.txt", sep = "\t")
#' ebseq("Row_count_3conditions.txt")
#'
#' @references Ning Leng and Christina Kendziorski (2020). EBSeq: An R package for gene and isoform differential expression analysis of RNA-seq data.
#' @param Row_count_matrix Row Count matrix txt file (Not normalized count matrix)
#' @export
#'
ebseq <- function(Row_count_matrix){
  DataMat <- data.matrix(read.table(Row_count_matrix, header = T, row.names = 1))
  collist <- gsub("\\_.+$", "", colnames(DataMat))
  if (length(unique(collist)) == 2) {
    name <- paste0(paste0(unique(collist)[1], "-vs-"), paste0(unique(collist)[2], "_EBseq"))}
  if (length(unique(collist)) == 3) {
    name <- paste0(paste0(unique(collist)[1], "-vs-"),
                   paste0(paste0(unique(collist)[2], "-vs-"),
                          paste0(unique(collist)[3], "_EBseq")))}

  print(name)
  vec <- c()
  for (i in 1:length(unique(collist))) {
    num <- length(collist[collist == unique(collist)[i]])
    vec <- c(vec, num)
  }
  ngvector <- NULL
  conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
  Sizes <- MedianNorm(DataMat)
  NormMat <- GetNormalizedMat(DataMat, Sizes)
  output_file <-paste0("result_of_", name)
  if (length(unique(collist)) == 2) {
    EBOut <- NULL
    EBOut <- EBTest(Data = DataMat, NgVector = ngvector, Conditions = conditions, sizeFactors = Sizes, maxround = 5)
    stopifnot(!is.null(EBOut))

    PP <- as.data.frame(GetPPMat(EBOut))
    fc_res <- PostFC(EBOut)

    results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,unlist(EBOut[["Mean"]][,1])[rownames(PP)], unlist(EBOut[["Mean"]][,2])[rownames(PP)])
    colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean")
    results <- results[order(results[,"PPDE"], decreasing = TRUE),]
    write.table(results, file = paste0(output_file, ".txt"), sep = "\t")

  } else {
    patterns <- GetPatterns(conditions)
    eename <- rownames(patterns)[which(rowSums(patterns) == length(unique(collist)))]
    stopifnot(length(eename) == 1)

    MultiOut <- NULL
    MultiOut <- EBMultiTest(Data = DataMat, NgVector = ngvector, Conditions = conditions, AllParti = patterns, sizeFactors = Sizes, maxround = 5)
    stopifnot(!is.null(MultiOut))

    MultiPP <- GetMultiPP(MultiOut)

    PP <- as.data.frame(MultiPP$PP)
    pos <- which(names(PP) == eename)
    probs <- rowSums(PP[,-pos])

    results <- cbind(PP, MultiPP$MAP[rownames(PP)], probs)
    colnames(results) <- c(colnames(PP), "MAP", "PPDE")
    ord <- order(results[,"PPDE"], decreasing = TRUE)
    results <- results[ord,]
    write.table(results, file = paste0(output_file, ".txt"), sep = "\t")

    write.table(MultiPP$Patterns, file = paste(output_file, ".pattern", sep = ""), sep = "\t")

    MultiFC <- GetMultiFC(MultiOut)
    write.table(MultiFC$CondMeans[ord,], file = paste(output_file, ".condmeans", sep = ""), sep = "\t")
  }
  norm_out_file <- paste0("Normalized_count_matrix_from_", paste0(name, ".txt"))
  write.table(NormMat, file = norm_out_file, sep = "\t")
}

