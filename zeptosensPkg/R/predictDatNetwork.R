#' predict data-driven only network.
#'
#' @param data input expression data frame. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param cutOff Manually Setup cut off value for the strength of edge. Default at 0.1.
#' @param nProt Number of proteins included in the drug perturbation data.
#' @param proteomicResponses Drug Perturbation data for analysis.
#' @param maxDist TBA (default: 1)
#' @return Parameter of regulization decided lowest BIC.Including regularize parameter(L1 norm parameter) as rho.
#' @examples
#' optimizeParameter(data = GeneExpresssion, prior = Priorindormation)
#' @concept zeptosensPkg
#' @export
predictDatNetwork <- function(data, cutOff = 0.1, nProt, proteomicResponses, maxDist = 1) {
  covmatrix <- cov(data)

  # optimize penalty parameter rho
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

  # Estimated inverse covariance (precision matrix)
  tmp <- glasso(covmatrix, rho = parameter)
  sigmaMatrix <- tmp$wi
  niter <- tmp$niter
  print(niter) # if niter = 10,000
  if (niter == 10000) {
    stop("ERROR: Algorithmn does not convergence!")
  }
  pcorMatrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      pcorMatrix[i, j] <- -sigmaMatrix[i, j] / sqrt(sigmaMatrix[i, i] * sigmaMatrix[j, j])
    }
  }

  t.edges <- pcorMatrix

  # cut off small edge value edges
  t.net <- as.data.frame(ifelse(abs(t.edges) >= cutOff & row(t.edges) != col(t.edges), t.edges, 0))
  colnames(t.net) <- colnames(data)
  rownames(t.net) <- colnames(data)

  # sum of edges
  nedges <- sum(t.net != 0)

  # Network to edgelist
  edgelist <- zeptosensPkg:::createSifFromMatrix(t.net = t.net, genelist = colnames(t.net))

  # network2 function into networkinfor
  networkInferred <- zeptosensPkg:::network2(
    wk = t.net, nProt = nProt,
    proteomicResponses = proteomicResponses,
    maxDist = maxDist
  )
  wk <- networkInferred$wk
  wks <- networkInferred$wks
  distInd <- networkInferred$distInd
  inter <- networkInferred$inter
  result <- list(
    rho = rho, nedges = nedges, t.net = t.net, edgelist = edgelist, bic = bic,
    wk = wk, wks = wks, distInd = distInd, inter = inter
  )
  return(result)
}
