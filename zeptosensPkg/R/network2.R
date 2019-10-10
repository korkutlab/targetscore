#' Network for adjusted-glasso
#' @param wk The upload network constructed In matrix form.
#' @param nProt number of proteins tested.
#' @param proteomicResponses TBA
#' @param distFile A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#'
#' @examples
#'
#' @concept zeptosensPkg
#'
network2 <- function(wk, nProt, proteomicResponses, maxDist, antibodyMapFile = NULL, distFile = NULL, verbose = F) {
  if (nProt != ncol(proteomicResponses)) {
    stop("ERROR: nProt is not equal to proteomicResponses column number")
  }

  # distInd(upstream,downstream)
  distInd <- matrix(Inf,
    ncol = nProt, nrow = nProt,
    dimnames = list(colnames(wk), colnames(wk))
  )
  for (i in 1:nProt) {
    for (j in 1:nProt) {
      if (wk[i, j] != 0) {
        distInd[i, j] <- 1
      } else {
        distInd[i, j] <- Inf
      }
    }
  }


  write.table(distInd, file = "distInd.txt", quote = F)
  wk <- (wk / max(abs(wk))) * maxDist
  wks <- wk

  inter <- (which(wk != 0, arr.ind = T))
  print(inter)
  networkInferred <- list(wk = wk, wks = wks, distInd = distInd, inter = inter)

  return(networkInferred)
}
