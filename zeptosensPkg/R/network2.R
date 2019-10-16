#' Network for adjusted-glasso
#'
#' @param wk Inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate and -1 for down regulate.
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#' @param max_dist Maximum edge strength value.(Default at 1)
#' @param antibody_map_file A list of antibodies, their associated genes, modification sites and effect.
#' @param dist_file A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#' @param verbose logical, whether to show additional debugging information
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

  wk <- (wk / max(abs(wk))) * max_dist
  wks <- wk

  inter <- (which(wk != 0, arr.ind = TRUE))
  print(inter)
  network_inferred <- list(wk = wk, wks = wks, dist_ind = dist_ind, inter = inter)

  return(network_inferred)
}
