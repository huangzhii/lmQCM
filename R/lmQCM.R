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
#' @param normalization Determine if normalization is needed on massive correlation coefficient matrix.
#' @return mergedCluster - An merged clusters group
#'
#' @examples
#' library(lmQCM)
#' library(Biobase)
#' data(sample.ExpressionSet)
#' data = assayData(sample.ExpressionSet)$exprs
#' lmQCM(data)
#'
#' @import genefilter
#' @import Biobase
#' @import nnet
#' @import stats
#' @export

lmQCM <- function(data_in,gamma=0.55,t=1,lambda=1,beta=0.4,minClusterSize=10,CCmethod="pearson", normalization = F) {
  message("Calculating massive correlation coefficient ...")
  cMatrix <- cor(t(data_in), method = CCmethod)
  diag(cMatrix) <- 0

  if(normalization){
    # Normalization
    cMatrix <- abs(cMatrix)
    D <- rowSums(cMatrix)
    D.half <- 1/sqrt(D)

    cMatrix <- apply(cMatrix, 2, function(x) x*D.half )
    cMatrix <- t(apply(cMatrix, 1, function(x) x*D.half ))
  }

  C <- localMaximumQCM(cMatrix, gamma, t, lambda)
  mergedCluster <- merging_lmQCM(C, beta, minClusterSize)
  message("Done.")
  return(mergedCluster)
}
