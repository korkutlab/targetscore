#' predict data-driven only network.
#'
#' @param data  input proteomics dataset for network inference. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param cut_off Manually Setup cut off value for the strength of edge. Default at 0.1.
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#'
#' @return Parameter of regulization decided lowest BIC.Including regularize parameter(L1 norm parameter) as rho.
#'
#' @importFrom glasso glasso
#'
#' @concept zeptosensPkg
#' @export
predict_dat_network <- function(data, cut_off = 0.1, n_prot, proteomic_responses) {
  covmatrix <- cov(data)

  # optimize penalty parameter rho
  rho <- seq(0.01, 1, length = 100)
  bic <- rho
  g_result <- NULL
  p_off_d <- NULL
  for (i in 1:100) {
    g_result <- glasso::glasso(covmatrix, rho[i])
    p_off_d <- sum(g_result$wi != 0 & col(covmatrix) < row(covmatrix))
    bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(data))
  }
  parameter <- rho[which.min(bic)]

  # Estimated inverse covariance (precision matrix)
  tmp <- glasso::glasso(covmatrix, rho = parameter)
  sigma_matrix <- tmp$wi
  niter <- tmp$niter
  print(niter) # if niter = 10,000
  if (niter == 10000) {
    stop("ERROR: Algorithmn does not convergence!")
  }
  pcor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      pcor_matrix[i, j] <- -sigma_matrix[i, j] / sqrt(sigma_matrix[i, i] * sigma_matrix[j, j])
    }
  }

  t_edges <- pcor_matrix

  # cut off small edge value edges
  t_net <- as.data.frame(ifelse(abs(t_edges) >= cut_off & row(t_edges) != col(t_edges), t_edges, 0))
  colnames(t_net) <- colnames(data)
  rownames(t_net) <- colnames(data)

  # sum of edges
  nedges <- sum(t_net != 0)

  # Network to edgelist
  edgelist <- zeptosensPkg::create_sif_from_matrix(t_net = t_net, genelist = colnames(t_net))

  # network2 function into networkinfor
  network_inferred <- zeptosensPkg::network2(
    wk = t_net, n_prot = n_prot,
    proteomic_responses = proteomic_responses,
  )
  wk <- network_inferred$wk
  wks <- network_inferred$wks
  dist_ind <- network_inferred$dist_ind
  inter <- network_inferred$inter

  result <- list(
    rho = rho, nedges = nedges, t_net = t_net, edgelist = edgelist, bic = bic,
    wk = wk, wks = wks, dist_ind = dist_ind, inter = inter
  )

  return(result)
}
