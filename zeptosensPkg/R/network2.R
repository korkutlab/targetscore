#' Network for adjusted-glasso
#'
#' @param wk Inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate and -1 for down regulate.
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#' @param max_dist Maximum edge strength value.(Default at 1)
#' @param dist_file A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#' @param verbose logical, whether to show additional debugging information
#'
#' @importFrom utils write.table
#'
#' @concept zeptosensPkg
#' @export
network2 <- function(wk, n_prot, proteomic_responses, max_dist, dist_file = NULL, verbose = FALSE) {
  if (n_prot != ncol(proteomic_responses)) {
    stop("ERROR: n_prot is not equal to proteomic_responses column number")
  }

  # converting the network to proteomic responces seq
  network <- wk
  protein_net <- data.frame(matrix(0, nrow = n_prot, ncol = n_prot))
  colnames(protein_net) <- colnames(proteomic_responses)
  rownames(protein_net) <- colnames(proteomic_responses)
  index <- which(colnames(proteomic_responses) %in% colnames(network))
  protein_net[index, index] <- network

  wk <- protein_net
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
