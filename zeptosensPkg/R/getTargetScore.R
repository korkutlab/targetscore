#' Target Score Pilot Run
#'
#' @param wk TBA
#' @param wks TBA
#' @param distInd TBA
#' @param inter TBA
#' @param nDose TBA
#' @param nProt TBA
#' @param proteomicResponses TBA
#' @param maxDist TBA (default: 1)
#' @param cellLine TBA
#' @param targetScoreOutputFile a filename to write total target score results (default: NULL)
#' @param matrixWkOutputFile TBA
#' @param targetScoreQValueFile a file namme to write statistical significance levels (default: NULL)
#' @param targetScorePValueFile a file namme to write statistical significance levels (default: NULL)
#' @param targetScoreDoseFile a filename to write dose dependent target score results (default: NULL)
#' @param nPerm number of random TS calculations for building the null distribution
#' @param verbose a flag for debugging output
#' @param tsFactor a scaling factor for the pathway component of the target score
#' @param fsFile a file with the functional score data
#' @param signedMatrixWkOutputFile TBA
#'
#' @details
#' data: multiple dose single drug perturbation
#' ts: integral_dose(fs*(xi+sigma_j(2^p*xj*product_k(wk))))
#' to be converted to integral_dose(fs*((xi/(stdev_xi)+sigma_j(2^p*((xj/stdev_xj)*product_k(wk))))
#' missing: For phosp and dephosp based wk, there is no 'exact match' between known and measured phospho-sites
#'
#' @examples
#'
#' @concept zeptosensPkg
#' @export
getTargetScore <- function(wk, wks, distInd, inter, nDose, nProt, proteomicResponses, 
                           maxDist = 1, nPerm, cellLine, targetScoreOutputFile = NULL,
                           matrixWkOutputFile = NULL, targetScoreQValueFile = NULL, 
                           targetScoreDoseFile = NULL, randomTargetScoreFile = NULL, 
                           targetScorePValueFile = NULL, verbose = TRUE, tsFactor = 1, 
                           fsFile, signedMatrixWkOutputFile = NULL) {

  # CALCULATE TARGET SCORE ----
  results <- calcTargetScore(
    wk = wk,
    wks = wks,
    distInd = distInd,
    inter = inter,
    nDose = nDose,
    nProt = nProt,
    proteomicResponses = proteomicResponses,
    maxDist = maxDist,
    cellLine = cellLine,
    verbose = TRUE,
    tsFactor = tsFactor,
    fsFile = fsFile
  )
  ts <- results$ts
  wk <- results$wk
  tsd <- results$tsd
  wks <- results$wks

  # random TS for each node over n permutations comes from randTargetScore.R
  randTs <- matrix(0, nrow = nProt, ncol = nPerm) 
  
  # p value for a given target score computed over the distribution from randTS
  pts <- matrix(0, ncol = 1, nrow = nProt) 

  # CREATE Q-VALUES ----
  for (k in 1:nPerm) {
    #        if(verbose) {
    cat("Permutation Iteration: ", k, "\n")
    #        }

    # print(fs) randomize the readouts over proteomic entities
    randProteomicResponses <- proteomicResponses

    # for(j in 1:ncol(randProteomicResponses))randProteomicResponses[,j] <- sample
    # (proteomicResponses[,j])
    for (i in 1:nrow(randProteomicResponses)) randProteomicResponses[i, ] <- sample(proteomicResponses[i, ])

    randTs[, k] <- calcTargetScore(
      wk = wk,
      wks = wks,
      distInd = distInd,
      inter = inter,
      nDose = nDose,
      nProt = nProt,
      proteomicResponses = randProteomicResponses,
      maxDist = maxDist,
      cellLine = cellLine,
      verbose = verbose,
      tsFactor = tsFactor,
      fsFile = fsFile
    )$ts

    # randTs[,k] <- as.matrix(rants) print('resi') print(resi$ts) randTs[,k]
  }

  for (i in 1:nProt) {
    mean <- mean(randTs[i, 1:nPerm])
    stdev <- sd(randTs[i, 1:nPerm])
    zval <- (ts[i] - mean) / (stdev)
    pts[i] <- 2 * pnorm(-abs(zval)) # pnorm(ts[i], mean = mean(randTs[i, 1:nPerm]), sd = sd(randTs[i, 1:nPerm]))

    if (verbose) {
      print(pts[i])
    }
  }
  q <- as.matrix(p.adjust(pts, method = "fdr", n = nProt))
  rownames(q) <- colnames(proteomicResponses)
  colnames(q) <- "FDR_adjusted_p"

  # WRITE OUTPUTS ----
  if (!is.null(matrixWkOutputFile)) {
    write.table(wk, file = matrixWkOutputFile, quote = FALSE, sep = "\t")
  }
  if (!is.null(signedMatrixWkOutputFile)) {
    write.table(wks, file = signedMatrixWkOutputFile, quote = FALSE, sep = "\t")
  }
  if (!is.null(targetScoreOutputFile)) {
    write.table(ts, file = targetScoreOutputFile, quote = FALSE, col.names = FALSE, sep = "\t")
  }

  if (!is.null(targetScoreDoseFile)) {
    write.table(tsd, file = targetScoreDoseFile, quote = FALSE, sep = "\t")
  }

  if (!is.null(randomTargetScoreFile)) {
    write.table(data.frame(randTs), file = randomTargetScoreFile, quote = FALSE, sep = "\t")
  }

  if (!is.null(targetScoreQValueFile)) {
    write.table(q, file = targetScoreQValueFile, quote = FALSE, sep = "\t")
  }
  rownames(randTs) <- colnames(proteomicResponses)
  rownames(pts) <- colnames(proteomicResponses)

  write.table(randTs, file = "randts.txt", quote = FALSE, sep = "\t")
  write.table(pts, file = targetScorePValueFile, quote = FALSE, sep = "\t")

  # RETURN RESULTS ----
  results <- list(ts = ts, wk = wk, tsd = tsd, q = q, wks = wks, pts = pts)

  return(results)
}
