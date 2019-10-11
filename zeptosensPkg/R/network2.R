#' Network for adjusted-glasso
#' @param wk The upload network constructed In matrix form.
#' @param n_prot number of proteins tested.
#' @param proteomic_responses TODO
#' @param dist_file A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#'
#' @importFrom utils write.table
#'
#' @concept zeptosensPkg
#' @export
network2 <- function(wk, n_prot, proteomic_responses, max_dist, antibody_map_file = NULL, 
                     dist_file = NULL, verbose = FALSE) {
  
  if (n_prot != ncol(proteomic_responses)) {
    stop("ERROR: n_prot is not equal to proteomic_responses column number")
  }

  # dist_ind(upstream,downstream)
  dist_ind <- matrix(Inf,
    ncol = n_prot, nrow = n_prot,
    dimnames = list(colnames(wk), colnames(wk))
  )
  for (i in 1:n_prot) {
    for (j in 1:n_prot) {
      if (wk[i, j] != 0) {
        dist_ind[i, j] <- 1
      } else {
        dist_ind[i, j] <- Inf
      }
    }
  }

  write.table(dist_ind, file = "dist_ind.txt", quote = FALSE)
  wk <- (wk / max(abs(wk))) * max_dist
  wks <- wk

  inter <- (which(wk != 0, arr.ind = TRUE))
  print(inter)
  network_inferred <- list(wk = wk, wks = wks, dist_ind = dist_ind, inter = inter)

  return(network_inferred)
}
