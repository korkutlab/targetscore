#' Compute a TargetScore for randomized data/compute P value for each TS given a network topology
#'
#' @param wk Inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate and -1 for down regulate.
#' @param wks Inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate, -1 for down regulate,
#' 2 for phosphorylation and -2 for dephosphorylation.
#' @param dist_ind A distance file of edgelist with a third column as the network distance
#'   between the genes in the interaction
#' @param edgelist Edgelist of inferred network
#' @param n_dose Dose number of input data.
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#' @param verbose a boolean to show debugging information
#' @param fs_dat Functional score file. A tab-delmited file with a header, each row is an
#'   antibody in the first column and functional score in the second column
#'   (i.e. 1 oncogene, 0 tumor supressor/oncogene, -1 tumor supressor characteristics)
#' @param dist_file A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#' @param ts_pathway_scale a scaling factor for the pathway component in the TargetScore
#'
#' @return a list is returned with the following entries:
#' {ts}{TargetScore values summed over individual drug doses}
#' {tsd}{TargetScore values for individual drug doses}
#' {wk}{inferred network matrix form with edge strength value estimated as the partial correlation.}
#' {wks}{inferred network matrix form with edge strength value estimated as the partial correlation.}
#'
#' @details
#' data: multiple dose single drug perturbation
#' ts: integral_dose(fs*(xi+sigma_j(2^p*xj*product_k(wk))))
#' missing: For phosp and dephosp based wk, there is no 'exact match' between known and measured phospho-sites
#'
#' @concept targetscore
#' @export
calc_target_score <- function(wk, wks, dist_ind, edgelist, n_dose, n_prot, proteomic_responses, fs_dat,
                              verbose = TRUE, ts_pathway_scale = 1, dist_file = NULL) {
  if(verbose) {
    tmp <- paste(capture.output(head(fs_dat, 3)), collapse = "\n")
    message("MSG: Functional score data (head):\n", tmp, "\n")
  }
  
  fs <- fs_dat

  # Calculate TS for each dose
  # print(n_dose)
  ## This will be the targetscore for each dose
  tsd <- matrix(0,
    nrow = n_dose, 
    ncol = n_prot,
    dimnames = list(
      rownames(proteomic_responses),
      colnames(proteomic_responses)
    )
  )
  
  ## This will be the pathway component for each targetscore
  tsp <- array(0:0,
    dim = c(n_dose, n_prot, n_prot),
    dimnames = list(
      rownames(proteomic_responses),
      colnames(proteomic_responses), 
      colnames(proteomic_responses)
    )
  )
  
  ## This will be the final targetscore summed over multiple doses; always 1 row
  ts <- matrix(0,
    ncol = n_prot, 
    nrow = 1,
    dimnames = list("targetScore", colnames(proteomic_responses))
  )

  edges_used <- 0
  
  for (i in 1:n_dose) {
    # downstream (target)
    for (j in 1:nrow(edgelist)) {
      # upstream
      # for (k in 1:n_prot) {
      node1 <- edgelist$source_node[j]
      node2 <- edgelist$target_node[j]
      # tsp[i,k,j] <- ts_pathway_scale*(2^-(dist_ind[k, j])) * proteomic_responses[i, k] * wk[k, j]
      
      if(node1 %in% colnames(proteomic_responses) & node2 %in% colnames(proteomic_responses)) {
        tsp[i, node1, node2] <- ts_pathway_scale * (2^-(dist_ind[node1, node2])) * proteomic_responses[i, node1] * wk[node1, node2]
        
        edges_used <- edges_used + 1
      }
    }
  }
  
  if(verbose) {
    message("MSG: Edges used per dose: ", floor(edges_used / n_dose), "\n")  
    
    network_nodes <- unique(c(edgelist$source_node, edgelist$target_node))
    response_nodes <- colnames(proteomic_responses)
    
    tmp <- paste(sort(setdiff(network_nodes, response_nodes)), collapse="|")
    message("MSG: IDs not present in network: ", tmp, "\n")
  }
  
  if(verbose) {
    debug <- data.frame(i=numeric(0), j=numeric(0), fs=numeric(0), self=numeric(0), network=numeric(0))
  } else {
    debug <- NULL
  }
  
  for (i in 1:n_dose) {
    for (j in 1:n_prot) {
      tsd[i, j] <- fs[j, 2] * (proteomic_responses[i, j] + sum(tsp[i, 1:n_prot, j]))
      
      if(verbose) {
        fs_tmp <- fs[j, 2]
        self_tmp <- proteomic_responses[i, j]
        network_tmp <- sum(tsp[i, 1:n_prot, j])
        
        tmp <- data.frame(i=i, j=j, fs=fs_tmp, self=self_tmp, network=network_tmp)
        debug <- rbind(debug, tmp)
      }
    }
  }

  ts <- colSums(tsd)
  # colnames(ts) <- colnames(proteomic_responses) rownames(ts) <- rownames(proteomic_responses)
  
  results <- list(wk = wk, wks = wks, ts = ts, tsd = tsd, debug=debug)
  return(results)
}
