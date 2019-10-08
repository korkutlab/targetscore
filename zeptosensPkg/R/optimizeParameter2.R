#' Choose the optimal regulization parameter for only data-driven network.
#'
#' @param data input expression data frame. Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @return Parameter of regulization decided lowest BIC.Including regularize parameter(L1 norm parameter) as rho.
#' @examples
#' optimizeParameter(data = GeneExpresssion, prior = Priorindormation)
#' @concept zeptosensPkg
#' @export
optimizeParameter2 <- function(data) {
  Covmatrix <- cov(data)
  rho <- seq(0.01, 1, length = 100)
  bic <- rho
  g.result <- c()
  p_off_d <- c()
  for (i in 1:100) {
    g.result <- glasso(Covmatrix, rho[i])
    p_off_d <- sum(g.result$wi != 0 & col(Covmatrix) < row(Covmatrix))
    bic[i] <- -2 * (g.result$loglik) + p_off_d * log(nrow(data))
  }
  parameter <- rho[which.min(bic)]
  return(parameters)
}
