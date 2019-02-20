setClass("QCMObject", representation(clusters.id = "list", clusters.names = "list",
                                     eigengene.matrix = "data.frame"))
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
#' @return QCMObject - An S4 Class with lmQCM results
#'
#' @examples
#' library(lmQCM)
#' library(Biobase)
#' data(sample.ExpressionSet)
#' data = assayData(sample.ExpressionSet)$exprs
#' data = fastFilter(data, 0.2, 0.2)
#' lmQCM(data)
#'
#' @import genefilter
#' @import Biobase
#' @import nnet
#' @import stats
#' @import methods
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

  clusters <- merging_lmQCM(C, beta, minClusterSize)
  # map rownames to clusters
  clusters.names = list()
  for (i in 1:length(clusters)){
    mc = clusters[[i]]
    clusters.names[[i]] = rownames(data_in)[mc]
  }
  # calculate eigengene
  eigengene.matrix <- matrix(0, nrow = length(clusters), ncol = dim(data_in)[2]) # Clusters * Samples

  for (i in 1:(length(clusters.names))) {
    geneID <- as.matrix(clusters.names[[i]])
    X <- data_in[geneID,]
    mu <- rowMeans(X)
    stddev <- rowSds(as.matrix(X), na.rm=TRUE) # standard deviation with 1/(n-1)
    XNorm <- sweep(X,1,mu) # normalize X
    XNorm <- apply(XNorm, 2, function(x) x/stddev)
    SVD <- svd(XNorm, LINPACK = FALSE)
    eigengene.matrix[i,] <- t(SVD$v[,1])
  }
  eigengene.matrix = data.frame(eigengene.matrix)
  colnames(eigengene.matrix) = colnames(data_in)

  QCMObject <- methods::new("QCMObject", clusters.id = clusters, clusters.names = clusters.names,
                   eigengene.matrix = eigengene.matrix)

  message("Done.")
  return(QCMObject)
}
