#' Create sif from network matrix.
#'
#' @param t_net p*p Gaussian gene Network matrix estimated.
#' @param row_genelist Gene name list corresponding to the row of gene Network matrix
#' @param col_genelist Gene name list corresponding to the column of gene Network matrix.
#'
#' @return Edgelist description data.frame contain G(V,E) including Vertex and Edges
#' (1 positive/-1 negative) edgevalue (the strength of association) and vertices numbers from the graph
#' 
#' @examples 
#' network_org <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
#' package = "targetscore"
#' ))
#' edgelist_wk <- create_sif_from_matrix(
#' t_net = network_org$wk,
#' col_genelist = colnames(network_org$wk),
#' row_genelist = rownames(network_org$wk)
#' )
#'
#' @concept targetscore
#' @export
create_sif_from_matrix <- function(t_net,
                                   row_genelist = rownames(t_net),
                                   col_genelist = colnames(t_net)) {
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
        node1[a] <- row_genelist[i]
        node1_n[a] <- i
        node2[a] <- col_genelist[j]
        node2_n[a] <- j
        edges_value[a] <- t_net[i, j]
        edges[a] <- ifelse(t_net[i, j] > 0, 1, -1)
        a <- a + 1
      }
    }
  }
  edgelist <- data.frame(
    source_node = node1,
    edge_value = edges,
    target_node = node2,
    edges_value = edges_value,
    source_index = node1_n,
    target_index = node2_n,
    stringsAsFactors = FALSE
  )
  return(edgelist)
}
