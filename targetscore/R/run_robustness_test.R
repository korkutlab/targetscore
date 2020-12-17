#' Run robustness test for data through the chosen algorithm.
#'
#' @param data input proteomics expression dataset for network inference. Coloumns as the gene, rows as the sample.With
#' colnames as the gene tags, rownames as the sample tags.
#' @param prior Prior information data frame ,with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags. Can be inferred from predict_bio_network() or any network resources.
#' Default at NULL and will be inferred from SignedPC database.
#' @param proteomic_responses RPPA data tested for drug pertubation.
#' @param n_prot Antibody number of input data.Here should be the number of antibody exist both in network
#' inference dataset and proteomic_responses dataset.
#' @param max_dist maximum distance between two antibody. (Default at 1)
#' @param mab_to_genes A list of antibodies, their associated genes, modification sites and effect.
#' @param boot_time Bootstrap time mannually set.Default at 100.
#' @param cut_off Manually set up cut off value for strength of edge. (Default at 0.1)
#' @param algorithm Flexible toolbox implementing network estimating algorithms for robustness test.
#' (\code{data_driven},or \code{hybrid_driven}).
#'
#' @return result of validation score. for random network and the network predicted with the algorithm.
#'
#' @concept targetscore
#' @export
run_robustness_test <- function(data, prior = NULL, proteomic_responses,
                                n_prot, max_dist = 1, boot_time = 100,
                                cut_off = 0.1, mab_to_genes,
                                algorithm = c("data_driven", "hybrid_driven")) {
  # Match values
  algorithm <- match.arg(algorithm)

  # Match the data
  index <- colnames(proteomic_responses[, which(colnames(proteomic_responses) %in% colnames(data))])
  data <- data[, index]
  data[is.na(data)] <- 0
  proteomic_responses <- proteomic_responses[, index]
  prior <- prior[index, index]

  # Prior extraction
  if (algorithm == "hybrid_driven" & is.null(prior)) {
    prior <- targetscore::predict_bio_network(
      n_prot = dim(proteomic_responses)[2],
      proteomic_responses = proteomic_responses,
      max_dist = 1,
      mab_to_genes = mab_to_genes
    )$wk
  }
  # Initial bootstrap result set up
  score_50p_tp <- array(0, dim = c(boot_time, 1))
  score_60p_tp <- array(0, dim = c(boot_time, 1))
  score_70p_tp <- array(0, dim = c(boot_time, 1))
  score_80p_tp <- array(0, dim = c(boot_time, 1))
  score_90p_tp <- array(0, dim = c(boot_time, 1))

  n_edges <- array(0, dim = c(boot_time, 1))
  n_edges_50p <- array(0, dim = c(boot_time, 1))
  n_edges_60p <- array(0, dim = c(boot_time, 1))
  n_edges_70p <- array(0, dim = c(boot_time, 1))
  n_edges_80p <- array(0, dim = c(boot_time, 1))
  n_edges_90p <- array(0, dim = c(boot_time, 1))

  # Boostrap Loop start here
  for (r in seq_len(boot_time)) {

    # Random selection of data
    data_50p <- data[sample(seq_len(nrow(data)), 0.5 * nrow(data), replace = FALSE), ]
    data_60p <- data[sample(seq_len(nrow(data)), 0.6 * nrow(data), replace = FALSE), ]
    data_70p <- data[sample(seq_len(nrow(data)), 0.7 * nrow(data), replace = FALSE), ]
    data_80p <- data[sample(seq_len(nrow(data)), 0.8 * nrow(data), replace = FALSE), ]
    data_90p <- data[sample(seq_len(nrow(data)), 0.9 * nrow(data), replace = FALSE), ]


    # Network construction
    if (algorithm == "data_driven") {
      # Data-driven Network
      network <- targetscore::predict_dat_network(
        data = data, cut_off = cut_off, n_prot = n_prot,
        proteomic_responses = proteomic_responses
      )$wk
      network_50p <- targetscore::predict_dat_network(
        data = data_50p, cut_off = cut_off, n_prot = n_prot,
        proteomic_responses = proteomic_responses
      )$wk
      network_60p <- targetscore::predict_dat_network(
        data = data_60p, cut_off = cut_off, n_prot = n_prot,
        proteomic_responses = proteomic_responses
      )$wk
      network_70p <- targetscore::predict_dat_network(
        data = data_70p, cut_off = cut_off, n_prot = n_prot,
        proteomic_responses = proteomic_responses
      )$wk
      network_80p <- targetscore::predict_dat_network(
        data = data_80p, cut_off = cut_off, n_prot = n_prot,
        proteomic_responses = proteomic_responses
      )$wk
      network_90p <- targetscore::predict_dat_network(
        data = data_90p, cut_off = cut_off, n_prot = n_prot,
        proteomic_responses = proteomic_responses
      )$wk
    }


    # Hybrid-driven network
    if (algorithm == "hybrid_driven") {
      network <- targetscore::predict_hybrid_network(
        data = data, cut_off = cut_off, n_prot = n_prot, prior = prior,
        proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
      )$wk
      network_50p <- targetscore::predict_hybrid_network(
        data = data_50p, cut_off = cut_off, n_prot = n_prot, prior = prior,
        proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
      )$wk
      network_60p <- targetscore::predict_hybrid_network(
        data = data_60p, cut_off = cut_off, n_prot = n_prot, prior = prior,
        proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
      )$wk
      network_70p <- targetscore::predict_hybrid_network(
        data = data_70p, cut_off = cut_off, n_prot = n_prot, prior = prior,
        proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
      )$wk
      network_80p <- targetscore::predict_hybrid_network(
        data = data_80p, cut_off = cut_off, n_prot = n_prot, prior = prior,
        proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
      )$wk
      network_90p <- targetscore::predict_hybrid_network(
        data = data_90p, cut_off = cut_off, n_prot = n_prot, prior = prior,
        proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
      )$wk
    }
    # Transform edge into(1/0/-1)
    network <- ifelse(network > cut_off, 1, ifelse(network < (-cut_off), -1, 0))
    network_50p <- ifelse(network_50p > cut_off, 1, ifelse(network_50p < (-cut_off), -1, 0))
    network_60p <- ifelse(network_60p > cut_off, 1, ifelse(network_60p < (-cut_off), -1, 0))
    network_70p <- ifelse(network_70p > cut_off, 1, ifelse(network_70p < (-cut_off), -1, 0))
    network_80p <- ifelse(network_80p > cut_off, 1, ifelse(network_80p < (-cut_off), -1, 0))
    network_90p <- ifelse(network_90p > cut_off, 1, ifelse(network_90p < (-cut_off), -1, 0))

    # Set score start value
    score_50p <- 0
    score_60p <- 0
    score_70p <- 0
    score_80p <- 0
    score_90p <- 0

    # Edges Re-driven & score count
    for (t in seq_len(length(index))) {
      for (p in seq_len(length(index))) {
        if (network_50p[t, p] == network[t, p] & network[t, p] != 0) {
          score_50p <- score_50p + 1
        }
        if (network_60p[t, p] == network[t, p] & network[t, p] != 0) {
          score_60p <- score_60p + 1
        }
        if (network_70p[t, p] == network[t, p] & network[t, p] != 0) {
          score_70p <- score_70p + 1
        }
        if (network_80p[t, p] == network[t, p] & network[t, p] != 0) {
          score_80p <- score_80p + 1
        }
        if (network_90p[t, p] == network[t, p] & network[t, p] != 0) {
          score_90p <- score_90p + 1
        }
      }
    }

    # Score re-store

    score_50p_tp[r] <- score_50p
    score_60p_tp[r] <- score_60p
    score_70p_tp[r] <- score_70p
    score_80p_tp[r] <- score_80p
    score_90p_tp[r] <- score_90p

    # Number of edges count
    n_edges[r] <- sum(network != 0)
    n_edges_50p[r] <- sum(network_50p != 0)
    n_edges_60p[r] <- sum(network_60p != 0)
    n_edges_70p[r] <- sum(network_70p != 0)
    n_edges_80p[r] <- sum(network_80p != 0)
    n_edges_90p[r] <- sum(network_90p != 0)
  }
  # loop end here,write out the result here with no overlap

  # Score Quatification & Summarize
  score_50p_mean <- mean(score_50p_tp / n_edges)
  data_50p_score <- score_50p_tp / n_edges

  score_60p_mean <- mean(score_60p_tp / n_edges)
  data_60p_score <- score_60p_tp / n_edges

  score_70p_mean <- mean(score_70p_tp / n_edges)
  data_70p_score <- score_70p_tp / n_edges

  score_80p_mean <- mean(score_80p_tp / n_edges)
  data_80p_score <- score_80p_tp / n_edges

  score_90p_mean <- mean(score_90p_tp / n_edges)
  data_90p_score <- score_90p_tp / n_edges

  result <- list(
    score_50p_mean = score_50p_mean, data_50p_score = data_50p_score, n_edges_50p = n_edges_50p,
    score_60p_mean = score_60p_mean, data_60p_score = data_60p_score, n_edges_60p = n_edges_60p,
    score_70p_mean = score_70p_mean, data_70p_score = data_70p_score, n_edges_70p = n_edges_70p,
    score_80p_mean = score_80p_mean, data_80p_score = data_80p_score, n_edges_80p = n_edges_80p,
    score_90p_mean = score_90p_mean, data_90p_score = data_90p_score, n_edges_90p = n_edges_90p,
    n_edges = n_edges
  )

  return(result)
}
