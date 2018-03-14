#' lmQCM: Main Routine for Gene Co-expression Analysis
#'
#' Author: Zhi Huang
#' @param data_in real-valued expression matrix with rownames indicating gene ID or gene symbol
#' @param gamma gamma value (default = 0.55)
#' @param t t value (default = 1)
#' @param lambda lambda value (default = 1)
#' @param beta beta value (default = 0.4)
#' @param minClusterSize minimum length of cluster to retain (default = 10)
#' @param CCmethod Methods for correlation coefficient calculation (default = "pearson"). Users can also pick "spearman".
#' @return mergedCluster
#' @import genefilter
#' @import Biobase
#' @importFrom nnet
#' @importFrom cor
#' @importFrom stats
#' @export

lmQCM <- function(data_in,gamma=0.55,t=1,lambda=1,beta=0.4,minClusterSize=10,CCmethod="pearson") {
  print("Calculating massive correlation coefficient ...")
  cMatrix = cor(t(data_in), method = CCmethod)
  diag(cMatrix) <- 0
  C <- localMaximumQCM(cMatrix, gamma, t, lambda)
  mergedCluster <- merging_lmQCM(C, beta, minClusterSize)
  print("Done.")
  return(mergedCluster)
}
