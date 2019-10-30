#' Choose the optimal regulization parameter and scale paramter for adjusted-glasso algorithm network construction
#'
#' @param data  input proteomics dataset for network inference. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param prior Prior information matrix of gene interaction, with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags.Can be inferred from Public data source (for example:SignedPC).
#'
#' @return Parameter list of regulization parameter decided by the prior information and the algorithmn lowest BIC.
#' Including regularize parameter(L1 norm parameter) as rho, scale parameter
#' (decided how much prior information contribute) as kappa, and regulization matrix
#' for the expression data as rho_m and Model's BIC matrix for differnet regularization parameters.
#'
#' @importFrom glasso glasso
#' @importFrom stats cov
#'
#' @concept zeptosensPkg
#' @export
optimize_parameter_hybrid <- function(data, prior) {
  index <- colnames(prior[, which(colnames(prior) %in% colnames(data))]) # match the data
  data <- data[, index]
  prior1 <- prior[index, index]

  # prior information extration
  prior1 <- ifelse(prior1 != 0, 1, 0) # information matrix of prior
  prior2 <- prior1 # symmetrical prior information
  for (i in 1:nrow(prior1)) {
    for (j in 1:ncol(prior1)) {
      if (prior1[i, j] != 0) {
        prior2[i, j] <- prior1[i, j]
        prior2[j, i] <- prior1[i, j]
      }
    }
  }
  prior2 <- ifelse(prior2 != 0, 1, 0)

  # getting the best tuning parameter from BIC minimization
  covmatrix <- cov(data)
  rho <- seq(0.01, 1, length = 100)
  bic <- matrix(NA, 100, 100)
  kappa <- rho
  rho_m <- NULL
  g_result <- NULL
  u <- matrix(1, nrow(prior2), ncol(prior2))
  p_off_d <- NULL
  for (i in 1:100) {
    for (j in 1:i) {
      rho_m <- rho[i] * u - kappa[j] * prior2
      g_result <- glasso::glasso(covmatrix, rho_m)
      p_off_d <- sum(g_result$wi != 0 & col(covmatrix) < row(covmatrix))
      bic[i, j] <- -2 * (g_result$loglik) + p_off_d * log(nrow(data))
      bic <- as.data.frame(bic)
      rownames(bic) <- rho
      colnames(bic) <- kappa
    }
  }
  pos <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)
  rho <- rho[pos[1]]
  kappa <- kappa[pos[2]]
  rho_m <- rho * u - kappa * prior2
  parameters <- list(rho_m = rho_m, rho = rho, kappa = kappa, bic = bic)


  return(parameters)
}
