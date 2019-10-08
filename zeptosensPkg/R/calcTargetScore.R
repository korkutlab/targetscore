#' Compute Target score for randomized data/compute P value for each TS given a network topology
#'
#' @param wk TBA
#' @param wks TBA
#' @param dist_ind TBA
#' @param inter TBA
#' @param nDose TBA
#' @param nProt TBA
#' @param proteomicResponses TBA
#' @param maxDist TBA (default: 1)
#' @param cellLine TBA
#' @param targetScoreOutputFile a filename to write target score results (default: NULL)
#' @param matrixWkOutputFile TBA
#' @param verbose a boolean to show debugging information
#' @param fsFile Functional score file. A tab-delmited file with a header, each row is an
#'   antibody in the first column and functional score in the second column
#'   (i.e. 1 oncogene, 0 tumor supressor/oncogene, -1 tumor supressor characteristics)
#' @param antibodyMapFile a listing of antibodies, their associated genes, and modification sites
#' @param distFile A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#' @param tsFactor a scaling factor for the pathway component in the target score
#'
#' @details
#' data: multiple dose single drug perturbation
#' ts: integral_dose(fs*(xi+sigma_j(2^p*xj*product_k(wk))))
#' missing: For phosp and dephosp based wk, there is no 'exact match' between known and measured phospho-sites
#'
#' @examples
#'
#' @concept zeptosensPkg
#' @export
calcTargetScore <- function(wk, wks, dist_ind, inter, nDose, nProt, proteomicResponses, maxDist = 1, cellLine, verbose = TRUE,
                            tsFactor = 1, fsFile, antibodyMapFile = NULL, distFile = NULL) {
  # LOAD & RANDOMIZE INTERNAL DATA ---- read function score
  # if(is.null(fsFile)) {
  #     fsFile <- system.file("targetScoreData", "fs.txt", package = "zeptosensPkg")
  # }

  fs <- read.table(fsFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

  if (verbose) {
    print(fs)
  }

  # calculate TS for each dose
  #    print(nDose)
  tsd <- matrix(0, nrow = nDose, ncol = nProt, dimnames = list(rownames(proteomicResponses), colnames(proteomicResponses)))
  tsp <- array(0:0, dim = c(nDose, nProt, nProt), dimnames = list(rownames(proteomicResponses), colnames(proteomicResponses), colnames(proteomicResponses)))
  ts <- matrix(0, ncol = nProt, nrow = 1, dimnames = list("targetScore", colnames(proteomicResponses)))

  for (i in 1:nDose) {
    # downstream (target)
    for (j in 1:nrow(inter)) {
      # upstream
      #            for (k in 1:nProt) {
      k <- as.numeric(inter[j, 1])
      l <- as.numeric(inter[j, 2])
      #        tsp[i,k,j] <- tsFactor*(2^-(dist_ind[k, j])) * proteomicResponses[i, k] * wk[k, j]
      tsp[i, k, l] <- tsFactor * (2^-(dist_ind[k, l])) * proteomicResponses[i, k] * wk[k, l]

      #            }
      #            tsd[i, j] <- fs[j, 2] * (proteomicResponses[i, j] + (tsp[i, j]))
    }
  }
  for (i in 1:nDose) {
    for (j in 1:nProt) {
      tsd[i, j] <- fs[j, 2] * (proteomicResponses[i, j] + sum(tsp[i, 1:nProt, j]))
    }
  }

  ts <- colSums(tsd)
  # colnames(ts) <- colnames(proteomicResponses) rownames(ts) <- rownames(proteomicResponses)
  results <- list(ts = ts, wk = wk, tsd = tsd, wks = wks)
  return(results)
}
