#' Create sif from network matrix.
#'
#' @param tNet p*p Gaussian gene Network matrix estimated.
#' @param genelist gene name list corresponding to the Gaussian gene Network estimated.
#' @return edgelist description data.frame contain G(V,E) including Vertex and Edges
#' (1 positive/-1 negative) edgevalue and vertex number of Gaussian Graphical Model.
#' @concept zeptosensPkg
#' @export
createSifFromMatrix <- function(tNet, genelist) {
  index <- genelist
  edgeNumber <- sum(tNet != 0)
  node1 <- array(0, dim = c(edgeNumber, 1))
  node1N <- array(0, dim = c(edgeNumber, 1))
  node2 <- array(0, dim = c(edgeNumber, 1))
  node2N <- array(0, dim = c(edgeNumber, 1))
  edges <- array(0, dim = c(edgeNumber, 1))
  edgesValue <- array(0, dim = c(edgeNumber, 1))
  a <- 1
  for (i in 1:nrow(tNet)) {
    for (j in 1:ncol(tNet)) {
      if (tNet[i, j] != 0) {
        node1[a] <- index[i]
        node1N[a] <- i
        node2[a] <- index[j]
        node2N[a] <- j
        edgesValue[a] <- tNet[i, j]
        edges[a] <- ifelse(tNet[i, j] > 0, 1, -1)
        a <- a + 1
      }
    }
  }
  edgelist <- as.data.frame(cbind(node1, edges, node2, edgesValue, node1N, node2N))
  colnames(edgelist) <- c("node1", "edges", "node2", "edge.value", "node1N", "node2N")
  return(edgelist)
}
