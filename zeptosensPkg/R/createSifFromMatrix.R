#' Create sif from network matrix.
#'
#' @param t_net p*p Gaussian gene Network matrix estimated.
#' @param genelist gene name list corresponding to the Gaussian gene Network estimated.
#' 
#' @return edgelist description data.frame contain G(V,E) including Vertex and Edges
#' (1 positive/-1 negative) edgevalue and vertex number of Gaussian Graphical Model.
#' 
#' @concept zeptosensPkg
#' @export
createSifFromMatrix <- function(t_net, genelist) {
  index <- genelist
  edge_number <- sum(t_net != 0)
  node1 <- array(0, dim = c(edge_number, 1))
  node1_n <- array(0, dim = c(edge_number, 1))
  node2 <- array(0, dim = c(edge_number, 1))
  node2_n <- array(0, dim = c(edge_number, 1))
  edges <- array(0, dim = c(edge_number, 1))
  edges_value <- array(0, dim = c(edge_number, 1))
  a <- 1
  for (i in 1:nrow(t_net)) {
    for (j in 1:ncol(t_net)) {
      if (t_net[i, j] != 0) {
        node1[a] <- index[i]
        node1_n[a] <- i
        node2[a] <- index[j]
        node2_n[a] <- j
        edges_value[a] <- t_net[i, j]
        edges[a] <- ifelse(t_net[i, j] > 0, 1, -1)
        a <- a + 1
      }
    }
  }
  edgelist <- as.data.frame(cbind(node1, edges, node2, edges_value, node1_n, node2_n))
  colnames(edgelist) <- c("node1", "edges", "node2", "edge.value", "node1_n", "node2_n")
  return(edgelist)
}
