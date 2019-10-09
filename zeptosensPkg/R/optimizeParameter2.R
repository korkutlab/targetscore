#' Choose the optimal regulization parameter for only data-driven network.
#'
#' @param data input expression data frame. Gene in coloumns and samples in row. 
#' With colnames as gene tags and rownames as sample tags.
#' @return Parameter of regulization decided lowest BIC.Including regularize parameter 
#' (L1 norm parameter) as rho.
#' @examples
#' optimizeParameter(data = GeneExpresssion, prior = Priorindormation)
#' @concept zeptosensPkg
#' @export
optimizeParameter2 <- function(data) {
  covmatrix <- cov(data)
  rho <- seq(0.01, 1, length = 100)
  bic <- rho
  gResult <- NULL
  pOffD <- NULL
  for (i in 1:100) {
    gResult <- glasso(covmatrix, rho[i])
    pOffD <- sum(gResult$wi != 0 & col(covmatrix) < row(covmatrix))
    bic[i] <- -2 * (gResult$loglik) + pOffD * log(nrow(data))
  }
  parameter <- rho[which.min(bic)]
  return(parameters)
}
