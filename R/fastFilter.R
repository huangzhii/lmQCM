#' fastFilter: Subroutine for filtering expression matrix
#'
#' Author: Zhi Huang
#' @param RNA an expression matrix (rows: genes; columns: samples)
#' @param lowest_percentile_mean a float value range 0-1
#' @param lowest_percentile_variance a float value range 0-1
#' @param var.func specify variance function
#' @return An filtered expression matrix
#' @import genefilter
#' @import Biobase
#' @import stats
#' @export
fastFilter <- function (RNA, lowest_percentile_mean = 0.2, lowest_percentile_variance = 0.2, var.func = "var"){
  RNA = as.matrix(RNA)
  rowIQRs <- function(eSet) {
    numSamp <- ncol(eSet)
    lowQ <- rowQ(eSet, floor(0.25 * numSamp))
    upQ <- rowQ(eSet, ceiling(0.75 * numSamp))
    upQ - lowQ
  }
  varFilter <- function (eset, var.cutoff = 0.5, filterByQuantile = TRUE, var.func = var.func) {
    if (deparse(substitute(var.func)) == "IQR") {
      message("Using row-wise IQR for calculating the variances.")
      vars <- rowIQRs(eset)
    } else {
      message("Calculating the variances.")
      vars <- apply(eset, 1, var.func)
    }
    if (filterByQuantile) {
      if (0 < var.cutoff && var.cutoff < 1) {
        quant = quantile(vars, probs = var.cutoff)
        selected = !is.na(vars) & vars > quant
      } else stop("Cutoff Quantile has to be between 0 and 1.")
    } else {
      selected <- !is.na(vars) & vars > var.cutoff
    }
    return(selected)
  }

  # Remove data with lowest m% mean exp value shared by all samples
  message("Note: For RNA data, we suppose input matrix (data frame) is with:")
  message("      Row: Genes;    Columns: Samples.")
  geneID = rownames(RNA)
  percentile = lowest_percentile_mean
  if (percentile > 0){
    RNAmean = apply(RNA, 1, mean)
    RNA_filtered1 = RNA[RNAmean > quantile(RNAmean, percentile), ]
    geneID_filtered1 = geneID[RNAmean > quantile(RNAmean, percentile)]
  } else {
    RNA_filtered1 = RNA
    geneID_filtered1 = geneID
  }
  message(sprintf("(%d genes, %d samples) after removing lowest %.2f%% mean expression value.",
                  dim(RNA_filtered1)[1], dim(RNA_filtered1)[2], percentile*100))

  # Remove data with lowest 10% variance across samples
  percentile = lowest_percentile_variance
  if (percentile > 0){
    if (dim(RNA_filtered1)[2] > 3){
      index <- varFilter(eset = RNA_filtered1, var.cutoff = percentile, var.func = var.func)
      RNA_filtered2 = RNA_filtered1[index, ]
      geneID_filtered2 = geneID_filtered1[index]
    } else{
      message("Cannot calculate order statistic on object with less than 3 columns, will not remove data based on variance.")
      RNA_filtered2 = RNA_filtered1
      geneID_filtered2 = geneID_filtered1
    }
  } else {
    RNA_filtered2 = RNA_filtered1
    geneID_filtered2 = geneID_filtered1
  }

  message(sprintf("(%d genes, %d samples) after removing lowest %.2f%% variance expression value.",
                  dim(RNA_filtered2)[1], dim(RNA_filtered2)[2], lowest_percentile_mean*100))
  return(RNA_filtered2)
}
