#' stdev over all samples
#' @param nX vertically merged data from all samples rows condition columns proteins
#' @param nSample number of samples
#' @param nProt number of proteins
#' @param nCond number of consitions
#' @concept zeptosensPkg
#' @export
#'
sampSdev <- function(nSample, nProt, nDose, nX) {



  #    nProt
  #    nSample
  #    nCondition
  sSdev <- array(0:0, dim = c(nDose, nProt))
  sSdev2 <- array(0:0, dim = c(nProt))

  #    for (i in 1:nSample) {

  for (j in 1:nProt) {
    for (k in 1:nDose) {
      sSdev[k, j] <- sd(nX[((k - 1) * nSample + 1):(k * nSample), j])
    }
  }
  #    }
  #    for (i in 1:nSample*nDose) {

  for (j in 1:nProt) {
    sSdev2[j] <- sd(nX[1:(nSample * nDose), j])
  }
  #    }
  return(sSdev2)
}
