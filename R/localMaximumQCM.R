#' localMaximumQCM: Subroutine for Creating Gene Clusters
#'
#' Author: Zhi Huang
#' @param cMatrix a correlation matirx
#' @param gamma gamma value (default = 0.55)
#' @param t t value (default = 1)
#' @param lambda lambda value (default = 1)
#' @return An unmerged clusters group 'C'
#' @import genefilter
#' @import Biobase
#' @import progress
#' @import stats
#' @export
localMaximumQCM <- function (cMatrix, gamma = 0.55, t = 1, lambda = 1){
  C <- list()
  nRow <- nrow(cMatrix)
  maxV <- apply(cMatrix, 2, max)
  maxInd <- apply(cMatrix, 2, which.max) # several diferrences comparing with Matlab results


  # Step 1 - find the local maximal edges
  # maxEdges <- matrix(0, nrow = 0, ncol = 2)
  # maxW <- matrix(0, nrow = 0, ncol = 1)
  # for (i in 1:nRow){
  #   if (maxV[i] == max(cMatrix[maxInd[i], ])) {
  #     maxEdges <- rbind(maxEdges, c(maxInd[i], i))
  #     maxW <- rbind(maxW, maxV[i])
  #   }
  # }
  lm.ind <- which(maxV == sapply(maxInd, function(x) max(cMatrix[x,])))
  maxEdges <- cbind(maxInd[lm.ind], lm.ind)
  maxW <- maxV[lm.ind]

  res <- sort.int(maxW, decreasing = TRUE, index.return=TRUE)
  sortMaxV <- res$x
  sortMaxInd <- res$ix
  sortMaxEdges <- maxEdges[sortMaxInd,]
  message(sprintf("Number of Maximum Edges: %d", length(sortMaxInd)))

  currentInit <- 1
  noNewInit <- 0

  nodesInCluster <- matrix(0, nrow = 0, ncol = 1)

  pb <- progress_bar$new(format = " Calculating [:bar] :percent eta: :eta",
                         total = length(sortMaxInd), clear = F, width=60)


  while ((currentInit <= length(sortMaxInd)) & (noNewInit == 0)) {
    pb$tick()
    if (sortMaxV[currentInit] < (gamma * sortMaxV[1]) ) {
      noNewInit <- 1
    }
    else {
      if ( (is.element(sortMaxEdges[currentInit, 1], nodesInCluster) == FALSE) & is.element(sortMaxEdges[currentInit, 2], nodesInCluster) == FALSE) {
        newCluster <- sortMaxEdges[currentInit, ]
        addingMode <- 1
        currentDensity <- sortMaxV[currentInit]
        nCp <- 2
        totalInd <- 1:nRow
        remainInd <- setdiff(totalInd, newCluster)
        # C = setdiff(A,B) for vectors A and B, returns the values in A that
        # are not in B with no repetitions. C will be sorted.
        while (addingMode == 1) {
          neighborWeights <- colSums(cMatrix[newCluster, remainInd])
          maxNeighborWeight <- max(neighborWeights)
          maxNeighborInd <- which.max(neighborWeights)
          c_v = maxNeighborWeight/nCp;
          alphaN = 1 - 1/(2*lambda*(nCp+t));
          if (c_v >= alphaN * currentDensity) {
            newCluster <- c(newCluster, remainInd[maxNeighborInd])
            nCp <- nCp + 1
            currentDensity <- (currentDensity*((nCp-1)*(nCp-2)/2)+maxNeighborWeight)/(nCp*(nCp-1)/2)
            remainInd <- setdiff(remainInd, remainInd[maxNeighborInd]);
          }
          else {
            addingMode <- 0
          }
        }
        nodesInCluster <- c(nodesInCluster, newCluster)
        C <- c(C, list(newCluster))
      }
    }
    currentInit <- currentInit + 1
  }
  message(" Calculation Finished.")
  return(C)
}
