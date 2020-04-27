#' Choose the optimal regulization parameter for only data-driven network.
#'
#' @param data input proteomics dataset for network inference. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param rho positive tuning parameter for elastic net penalty. Default at 10^seq(-2, 0, 0.02).
#'
#' @return Parameter of regulization BIC error calculated.
#' \item{rho}{penalty parameter optimized.}
#' \item{bic}{BIC error calculated for penalty parameters.}
#' 
#' @examples 
#' # read proteomic responce file
#' signaling_responses <- read.csv(system.file("test_data", "TCGA-BRCA-L4.csv", package = "zeptosensPkg"), row.names = 1)
#'  parameters <- zeptosensPkg::optimize_parameter_dat(data = signaling_responses)
#' 
#' @importFrom glasso glasso
#'
#' @concept zeptosensPkg
#' @export
optimize_parameter_dat <- function(data, rho = 10^seq(-2, 0, 0.02)) {

  # Calculate covariance matrix
  covmatrix <- (nrow(data) - 1) / nrow(data) * stats::cov(data)

  # Set initial value
  bic <- rho
  g_result <- NULL
  p_off_d <- NULL

  # Calculate BIC for grid-search
  for (i in seq_len(length(rho))) {
    g_result <- glasso::glasso(covmatrix, rho[i], nobs = nrow(covmatrix))
    p_off_d <- sum(g_result$wi != 0 & col(covmatrix) < row(covmatrix))
    bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(data))
  }

  rho <- rho[which.min(bic)]
  parameters <- list(rho = rho, bic = bic)
  return(parameters)
}
