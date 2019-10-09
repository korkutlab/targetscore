#' Predicted Directional GeneNetwork
#'
#' @param data input expression data. Coloumns as the gene, rows as the sample. 
#' With colnames as the gene tags, rownames as the sample tags.
#' @param prior prior information matrix with colnames and rownames as the gene tags. 
#' Can extracted from predictBioNetwork function or inferred from any resources.
#' @param rho regulization parameter
#' @param kappa scaler parameter
#' @return estimatedNetwork include list of estimated directional partial correlation 
#' gene network rho as the regulization parametwe and kappa as the scaler parameter.
#' @concept zeptosensPkg
#' @export
predictDirectionalNetwork <- function(data, prior, rho, kappa, cut.off) {
  index <- colnames(prior[, which(colnames(prior) %in% colnames(data))]) # match the data

  data <- data[, index]
  prior <- prior[index, index]
  prior <- ifelse(prior != 0, 1, 0) # information matrix of prior

  u <- matrix(1, ncol(data), ncol(data))
  rhoM <- rho * u - kappa * prior
  pc <- cov(data)
  # Network construction with directional prior information
  sigmaMatrix <- glasso(pc, rho = rhoM)$wi
  pcorMatrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      pcorMatrix[i, j] <- -sigmaMatrix[i, j] / sqrt(sigmaMatrix[i, i] * sigmaMatrix[j, j])
    }
  }

  tEdgesRhoadjusted <- pcorMatrix

  # get direction for the network as the bigger covariance estimated indicated the upper stream gene;
  tEdgesRhoadjustedD <- matrix(0, nrow = nrow(tEdgesRhoadjusted), ncol = ncol(tEdgesRhoadjusted))
  for (j in 1:ncol(tEdgesRhoadjusted)) {
    for (i in 1:nrow(tEdgesRhoadjusted)) {
      if (abs(tEdgesRhoadjusted[i, j]) > abs(tEdgesRhoadjusted[j, i])) {
        tEdgesRhoadjustedD[i, j] <- tEdgesRhoadjusted[i, j]
      }
      if (abs(tEdgesRhoadjusted[i, j]) < abs(tEdgesRhoadjusted[j, i])) {
        tEdgesRhoadjustedD[j, i] <- tEdgesRhoadjusted[j, i]
      }
      if (abs(tEdgesRhoadjusted[i, j]) == abs(tEdgesRhoadjusted[j, i])) {
        tEdgesRhoadjustedD[i, j] <- tEdgesRhoadjusted[i, j]
        tEdgesRhoadjustedD[j, i] <- tEdgesRhoadjusted[j, i]
      }
    }
  }
  # cutoff at cutoff point
  tEdgesRhoadjustedD <- ifelse(abs(tEdgesRhoadjustedD) > cut.off, tEdgesRhoadjustedD, 0)
  estimatedNetwork <- list(edge = tEdgesRhoadjustedD, rho, kappa)
  return(estimatedNetwork)
}
