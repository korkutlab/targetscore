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
#' @param percentage The percentage of data for training dataset from dataset source. Default at 50.
#' @param cut_off Manually set up cut off value for strength of edge. (Default at 0)
#'
#' @return result of validation score. for random network and the network predicted with the algorithm.
#'
#' @concept zeptosensPkg
#' @export
run_robustness_test <- function(data, prior = NULL, proteomic_responses, n_prot, max_dist = 1,
                                boot_time = 1000, cut_off = 0, fold = 5, mab_to_genes) {
  # Match the data
  index <- colnames(proteomic_responses[, which(colnames(proteomic_responses) %in% colnames(data))])
  data <- data[, index]
  data[is.na(data)] <- 0
  proteomic_responses <- proteomic_responses[, index]

  # Bootstrap result set up
  score_dat_tp <- array(0, dim = c(boot_time, 1))
  score_datr_tp <- array(0, dim = c(boot_time, 1))
  score_hyb_tp <- array(0, dim = c(boot_time, 1))
  score_hybr_tp <- array(0, dim = c(boot_time, 1))

  n_edges_dat <- array(0, dim = c(boot_time, 1))
  n_edges_datr <- array(0, dim = c(boot_time, 1))
  n_edges_hyb <- array(0, dim = c(boot_time, 1))
  n_edges_hybr <- array(0, dim = c(boot_time, 1))

  # Boostrap Loop start here
  for (r in 1:boot_time) {

    # Randomization of Data
    random_data <- as.data.frame(t(apply(data, 1, function(x) {
      sample(x, replace = FALSE)
    })))
    colnames(random_data) <- colnames(data)

    # Split the data and random data into training and valid

    valid_n <- sample(seq_len(nrow(data)), (1 / fold) * nrow(data), replace = FALSE) # default fold=5
    valid_data <- data[valid_n, ]
    train_data <- data[-valid_n, ]

    valid_nr <- sample(seq_len(nrow(random_data)), (1 / fold) * nrow(random_data), replace = FALSE) # default fold=5
    valid_datar <- data[valid_nr, ]
    train_datar <- data[-valid_nr, ]

    # Network construction

    # Data-driven Network
    # Training-data network
    train_network_dat <- zeptosensPkg::predict_dat_network(
      data = train_data, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses
    )$wk
    train_network_datr <- zeptosensPkg::predict_dat_network(
      data = train_datar, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses
    )$wk
    # Valid-data network
    valid_network_dat <- zeptosensPkg::predict_dat_network(
      data = valid_data, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses
    )$wk
    valid_network_datr <- zeptosensPkg::predict_dat_network(
      data = valid_datar, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses
    )$wk
    # Hybrid-driven network
    # Training-data network
    train_network_hyb <- zeptosensPkg::predict_hybrid_network(
      data = train_data, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
    )$wk
    train_network_hybr <- zeptosensPkg::predict_hybrid_network(
      data = train_datar, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
    )$wk
    # Valid-data network
    valid_network_hyb <- zeptosensPkg::predict_hybrid_network(
      data = valid_data, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
    )$wk
    valid_network_hybr <- zeptosensPkg::predict_hybrid_network(
      data = valid_datar, cut_off = cut_off, n_prot = n_prot,
      proteomic_responses = proteomic_responses, mab_to_genes = mab_to_genes
    )$wk

    # Transform edge into(1/0/-1)
    train_network_dat <- ifelse(train_network_dat > cut_off, 1, ifelse(train_network_dat < (-cut_off), -1, 0))
    train_network_datr <- ifelse(train_network_datr > cut_off, 1, ifelse(train_network_datr < (-cut_off), -1, 0))
    valid_network_dat <- ifelse(valid_network_dat > cut_off, 1, ifelse(valid_network_dat < (-cut_off), -1, 0))
    valid_network_datr <- ifelse(valid_network_datr > cut_off, 1, ifelse(valid_network_datr < (-cut_off), -1, 0))
    train_network_hyb <- ifelse(train_network_hyb > cut_off, 1, ifelse(train_network_hyb < (-cut_off), -1, 0))
    train_network_hybr <- ifelse(train_network_hybr > cut_off, 1, ifelse(train_network_hybr < (-cut_off), -1, 0))
    valid_network_hyb <- ifelse(valid_network_hyb > cut_off, 1, ifelse(valid_network_hyb < (-cut_off), -1, 0))
    valid_network_hybr <- ifelse(valid_network_hybr > cut_off, 1, ifelse(valid_network_hybr < (-cut_off), -1, 0))

    # Number of edges count
    n_edges_dat[r] <- sum(train_network_dat != 0)
    n_edges_datr[r] <- sum(train_network_datr != 0)
    n_edges_hyb[r] <- sum(train_network_hyb != 0)
    n_edges_hybr[r] <- sum(train_network_hybr != 0)

    # Set score start value
    score_dat <- 0
    score_datr <- 0
    score_hyb <- 0
    score_hybr <- 0

    # Edges Re-driven & score count
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (train_network_dat[t, p] == valid_network_dat[t, p] & train_network_dat[t, p] != 0) {
          score_dat <- score_dat + 1
        }
        if (train_network_hyb[t, p] == valid_network_hyb[t, p] & train_network_hyb[t, p] != 0) {
          score_hyb <- score_hyb + 1
        }
        if (train_network_datr[t, p] == valid_network_datr[t, p] & train_network_datr[t, p] != 0) {
          score_datr <- score_datr + 1
        }
        if (train_network_hybr[t, p] == valid_network_hybr[t, p] & train_network_hybr[t, p] != 0) {
          score_hybr <- score_hybr + 1
        }
      }
    }

    # Score re-store

    score_dat_tp[r] <- score_dat
    score_datr_tp[r] <- score_datr
    score_hyb_tp[r] <- score_hyb
    score_hybr_tp[r] <- score_hybr
  }
  # loop end here,write out the result here with no overlap

  # Score Quatification & Summarize
  # Data-driven for data
  score_dat_mean <- mean(score_dat_tp / n_edges_dat)
  data_net_score <- score_dat_tp / n_edges_dat

  # Data-driven for randomized data
  score_datr_mean <- mean(score_datr_tp / n_edges_datr)
  data_net_scorer <- score_datr / n_edges_datr

  # Hybrid-approach for data
  score_hyb_mean <- mean(score_hyb / n_edges_hyb)
  hyb_net_score <- score_hyb / n_edges_hyb

  # Hybrid-approach for randomized data
  score_hybr_mean <- mean(score_hybr / n_edges_hybr)
  hyb_net_scorer <- score_hybr / n_edges_hybr

  result <- list(
    score_dat_mean = score_dat_mean, data_net_score = data_net_score, n_edges_dat = n_edges_dat,
    score_datr_mean = score_datr_mean, data_net_scorer = data_net_scorer, n_edges_datr = n_edges_datr,
    score_hyb_mean = score_hyb_mean, hyb_net_score = hyb_net_score, n_edges_hyb = n_edges_hyb,
    score_hybr_mean = score_hybr_mean, hyb_net_scorer = hyb_net_scorer, n_edges_hybr = n_edges_hybr
  )

  return(result)
}
