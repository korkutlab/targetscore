#' Choose the optimal regulization parameter for only data-driven network.
#'
#' @param data input proteomics dataset for network inference. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#'
#' @return Parameter of regulization decided lowest BIC. Including regularize parameter
#' (L1 norm parameter) as rho.
#'
#' @importFrom glasso glasso
#'
#' @concept zeptosensPkg
#' @export
optimize_parameter_dat <- function(data) {
  covmatrix <- cov(data)
  rho <- seq(0.01, 1, length = 100)
  bic <- rho
  g_result <- NULL
  p_off_d <- NULL
  for (i in 1:100) {
    g_result <- glasso::glasso(covmatrix, rho[i])
    p_off_d <- sum(g_result$wi != 0 & col(covmatrix) < row(covmatrix))
    bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(data))
  }
  rho <- rho[which.min(bic)]
  parameters <- list(rho = rho, bic = bic)
  return(parameters)
}
