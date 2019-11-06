#' Choose the optimal regulization parameter for only data-driven network.
#'
#' @param data input proteomics dataset for network inference. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param rho positive tuning parameter for elastic net penalty. Default at seq(0.01,1,length=100).
#'
#' @return Parameter of regulization BIC error calculated.
#' \item{rho}{penalty parameter.}
#' \item{bic}{BIC error calculated for penalty parameters.}
#'
#' @importFrom glasso glasso
#'
#' @concept zeptosensPkg
#' @export
optimize_parameter_dat <- function(data, rho = seq(0.01, 1, length = 100)) {

  # Calculate covariance matrix
  covmatrix <- cov(data)

  # Set initial value
  bic <- rho
  g_result <- NULL
  p_off_d <- NULL

  # Calculate BIC for grid-search
  for (i in seq_len(length(rho))) {
    g_result <- glasso::glasso(covmatrix, rho[i])
    p_off_d <- sum(g_result$wi != 0 & col(covmatrix) < row(covmatrix))
    bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(data))
  }

  rho <- rho[which.min(bic)]
  parameters <- list(rho = rho, bic = bic)
  return(parameters)
}
