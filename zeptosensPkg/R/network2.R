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

  dist_ind <- matrix(Inf, ncol = nProt, nrow = nProt, dimnames = list(colnames(wk), colnames(wk))) # dist_ind(upstream,downstream)
  for (i in 1:nProt) {
    for (j in 1:nProt) {
      if (wk[i, j] != 0) {
        dist_ind[i, j] <- 1
      } else {
        dist_ind[i, j] <- Inf
      }
    }
  }


  write.table(dist_ind, file = "dist_ind.txt", quote = F)
  wk <- (wk / max(abs(wk))) * maxDist
  wks <- wk

  inter <- (which(wk != 0, arr.ind = T))
  print(inter)
  networkInferred <- list(wk = wk, wks = wks, dist_ind = dist_ind, inter = inter)

  return(networkInferred)
}
