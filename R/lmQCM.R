#' lmQCM: An Algorithm for Detecting Weak Quasi-Cliques in Weighted Graph with Applications in Gene Co-Expression Module Discovery in Cancers
#' 03/14/2018 1:59PM The Pi Day
#' Author: Zhi Huang
#' @param data_in, gamma, t, lambda, beta, minClusterSize, CCmethod
#' @return mergedCluster

library(genefilter)
library(Biobase)

localMaximumQCM <- function (cMatrix, gamma, t, lambda){
  C <- list()
  nRow <- nrow(cMatrix)
  maxV <- apply(cMatrix, 2, max)
  maxInd <- apply(cMatrix, 2, which.max) # several diferrences comparing with Matlab results
  # Step 1 - find the local maximal edges
  maxEdges <- matrix(0, nrow = 0, ncol = 2)
  maxW <- matrix(0, nrow = 0, ncol = 1)
  for (i in 1:nRow){
    if (maxV[i] == max(cMatrix[maxInd[i], ])) {
      maxEdges <- rbind(maxEdges, c(maxInd[i], i))
      maxW <- rbind(maxW, maxV[i])
    }
  }

  res <- sort.int(maxW, decreasing = TRUE, index.return=TRUE)
  sortMaxV <- res$x
  sortMaxInd <- res$ix
  sortMaxEdges <- maxEdges[sortMaxInd,]
  print(sprintf("Number of Maximum Edges: %d", length(sortMaxInd)))

  currentInit <- 1
  noNewInit <- 0

  nodesInCluster <- matrix(0, nrow = 0, ncol = 1)
  while ((currentInit <= length(sortMaxInd)) & (noNewInit == 0)) {
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
          maxNeighborInd <- which.is.max(neighborWeights)
          c_v = maxNeighborWeight/nCp;
          alphaN = 1 - 1/(2*lambda*(nCp+t));
          if (c_v >= alphaN * currentDensity) {
            newCluster <- c(newCluster, remainInd[maxNeighborInd])
            nCp <- nCp + 1
            currentDensity <- (currentDensity*((nCp-1)*(nCp-2)/2)+maxNeighborWeight)/(nCp*(nCp-1)/2)
            remainInd = setdiff(remainInd, remainInd[maxNeighborInd]);
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
  return(C)
}

merging_lmQCM <- function(C, beta, minClusterSize){
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

  print("start merge")
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
        print(sprintf("The length of merged Cluster: %d", length(mergedCluster)))
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
  return(mergedCluster)
}

lmQCM <- function(data_in,gamma,t,lambda,beta,minClusterSize,CCmethod="pearson") {
  print("Calculating massive correlation coefficient ...")
  cMatrix = cor(t(data_in), method = CCmethod)
  diag(cMatrix) <- 0
  C <- localMaximumQCM(cMatrix, gamma, t, lambda)
  mergedCluster <- merging_lmQCM(C, beta, minClusterSize)
  return(mergedCluster)
}
