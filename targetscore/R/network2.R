#' Network for adjusted-glasso
#'
#' @param wk Inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate and -1 for down regulate.
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#' @param dist_file A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#' @param verbose logical, whether to show additional debugging information
#'
#' @return a list is returned with the following entries:
#' * "wk" inferred network matrix form with edge strength value default at 1 for
#' upregulate and -1 for down regulate.
#' * "wks" as inferred network matrix form with edge strength value default at 1 for
#' upregulate and -1 for down regulate and 2 for phosphorylation and -2 for dephosphorylation.
#' * "dist_ind" A distance file of edgelist with a third column as the network distance
#' between the genes in the interaction.
#' * "inter" file as edgelist of inferred network.
#' 
#' @examples 
#' # Read proteomic response for cellline1
#' file <- system.file("test_data", "BT474.csv", package = "targetscore")
#' proteomic_responses <- read.csv(file, row.names = 1)
#'  
#'  # Read in network output
#'  wk_org <- readRDS(system.file("test_data_files", "predict_hybrid_network_network_output.rds",
#'  package = "targetscore"
#'  ))
#'  
#'  network <- targetscore::network2(
#'  wk <- wk_org$wk,
#'  n_prot = dim(proteomic_responses)[2],
#'  proteomic_responses = proteomic_responses
#'  )
#' 
#' @importFrom utils write.table
#'
#' @concept targetscore
#' @export
network2 <- function(wk, n_prot, proteomic_responses, dist_file = NULL, verbose = FALSE) {
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

  wk <- wk / max(abs(wk))
  wks <- wk

  inter <- (which(wk != 0, arr.ind = TRUE))
  print(inter)
  network_inferred <- list(wk = wk, wks = wks, dist_ind = dist_ind, inter = inter)

  return(network_inferred)
}
