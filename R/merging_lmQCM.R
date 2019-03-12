#' merging_lmQCM: Subroutine for Merging Gene Clusters
#'
#' Author: Zhi Huang
#' @param C Resulting clusters
#' @param beta beta value (default = 0.4)
#' @param minClusterSize minimum length of cluster to retain (default = 10)
#' @return mergedCluster - An merged clusters group
#' @import genefilter
#' @import Biobase
#' @import stats
#' @export

merging_lmQCM <- function(C, beta=0.4, minClusterSize=10){
  # Merge the overlapped networks
  sizeC <- matrix(0, nrow = 0, ncol = length(C))
  for (i in 1:length(C)){
    sizeC[i] <- length(C[[i]])
  }
  res <- sort.int(sizeC, decreasing = TRUE, index.return=TRUE)
  sortC <- res$x
  sortInd <- res$ix

  C <- C[sortInd] # Still C, but sorted based on number of elements in each cell

  ind <- which(sortC >= minClusterSize)

  mergedCluster <- C[ind]
  mergeOccur <- 1
  currentInd <- 0

  message(sprintf(" %d Modules before merging.", length(C)))
  while (mergeOccur == 1) {
    mergeOccur <- 0
    while (currentInd < length(mergedCluster)){
      currentInd <- currentInd + 1
      if (currentInd < length(mergedCluster)){
        keepInd <- 1:currentInd
        for (j in (currentInd+1) : length(mergedCluster)) {
          interCluster <- intersect(mergedCluster[[currentInd]], mergedCluster[[j]]);
          if (length(interCluster) >= beta*min(length(mergedCluster[[j]]), length(mergedCluster[[currentInd]]))) {
            mergedCluster[currentInd] <- list(union(mergedCluster[[currentInd]], mergedCluster[[j]]))
            mergeOccur <- 1
          }
          else {
            keepInd <- c(keepInd, j)
          }
        }
        mergedCluster <- mergedCluster[keepInd]
        # message(sprintf("The length of merged Cluster: %d", length(mergedCluster)))
      }
    }
    sizeMergedCluster <- matrix(0, nrow = 0, ncol = length(mergedCluster))
    for (i in 1 : length(mergedCluster)) {
      sizeMergedCluster[i] <- length(mergedCluster[[i]])
    }
    res <- sort.int(sizeMergedCluster, decreasing = TRUE, index.return=TRUE)
    sortSize <- res$x
    sortMergedInd <- res$ix
    mergedCluster <- mergedCluster[sortMergedInd]
    currentInd <- 0
  }
  for (i in 1:length(mergedCluster)){
    mergedCluster[[i]] <- unname(mergedCluster[[i]])
  }
  message(sprintf(" %d Modules remain after merging.", length(mergedCluster)))
  return(mergedCluster)
}
