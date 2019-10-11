#' Predicted Directional GeneNetwork
#'
#' @param data input expression data. Coloumns as the gene, rows as the sample.
#' With colnames as the gene tags, rownames as the sample tags.
#' @param prior prior information matrix with colnames and rownames as the gene tags.
#' Can extracted from predictBioNetwork function or inferred from any resources.
#' @param rho regulization parameter
#' @param kappa scaler parameter
#' @param cut_off TODO
#' 
#' @return estimated_network include list of estimated directional partial correlation
#' gene network rho as the regulization parametwe and kappa as the scaler parameter.
#' 
#' @concept zeptosensPkg
#' @export
predictDirectionalNetwork <- function(data, prior, rho, kappa, cut_off) {
  index <- colnames(prior[, which(colnames(prior) %in% colnames(data))]) # match the data

  data <- data[, index]
  prior <- prior[index, index]
  prior <- ifelse(prior != 0, 1, 0) # information matrix of prior

  u <- matrix(1, ncol(data), ncol(data))
  rho_m <- rho * u - kappa * prior
  pc <- cov(data)
  # Network construction with directional prior information
  sigma_matrix <- glasso(pc, rho = rho_m)$wi
  pcor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      pcor_matrix[i, j] <- -sigma_matrix[i, j] / sqrt(sigma_matrix[i, i] * sigma_matrix[j, j])
    }
  }

  t_edges_rhoadjusted <- pcor_matrix

  # get direction for the network as the bigger covariance estimated indicated the upper stream gene;
  t_edges_rhoadjusted_d <- matrix(0, nrow = nrow(t_edges_rhoadjusted), ncol = ncol(t_edges_rhoadjusted))
  for (j in 1:ncol(t_edges_rhoadjusted)) {
    for (i in 1:nrow(t_edges_rhoadjusted)) {
      if (abs(t_edges_rhoadjusted[i, j]) > abs(t_edges_rhoadjusted[j, i])) {
        t_edges_rhoadjusted_d[i, j] <- t_edges_rhoadjusted[i, j]
      }
      if (abs(t_edges_rhoadjusted[i, j]) < abs(t_edges_rhoadjusted[j, i])) {
        t_edges_rhoadjusted_d[j, i] <- t_edges_rhoadjusted[j, i]
      }
      if (abs(t_edges_rhoadjusted[i, j]) == abs(t_edges_rhoadjusted[j, i])) {
        t_edges_rhoadjusted_d[i, j] <- t_edges_rhoadjusted[i, j]
        t_edges_rhoadjusted_d[j, i] <- t_edges_rhoadjusted[j, i]
      }
    }
  }
  # cutoff at cutoff point
  t_edges_rhoadjusted_d <- ifelse(abs(t_edges_rhoadjusted_d) > cut_off, t_edges_rhoadjusted_d, 0)
  estimated_network <- list(edge = t_edges_rhoadjusted_d, rho, kappa)
  return(estimated_network)
}
