#' Run Cross validation for data through the chosen algorithm.
#'
#' @param data input expression data. Coloumns as the gene, rows as the sample.With
#' colnames as the gene tags, rownames as the sample tags.
#' @param prior Prior information data frame ,with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags. Can be inferred from predict_bio_network() or any network resources.
#' @param boot_time Bootstrap time mannually set.Default at 1000.
#' @param fold The fold for training and test dataset. Default at 5.
#' @param cut_off Manually set up cut off value for strength of edge. (Default at 0.1)
#'
#' @return result of validation score. for random network and the network predicted with the algorithm.
#'
#' @importFrom glasso glasso
#' @importFrom GeneNet ggm.simulate.pcor
#' @importFrom stats var
#'
#' @concept zeptosensPkg
#' @export
run_cross_validation <- function(data, prior, boot_time = 1000, cut_off = 0.1, fold = 5) {
  index <- colnames(prior[, which(colnames(prior) %in% colnames(data))]) # match the data

  data <- data[, index]
  data[is.na(data)] <- 0
  prior1 <- prior[index, index]
  prior1 <- ifelse(prior1 != 0, 1, 0) # turn prior information into binary data

  prior2 <- prior1 # symmetrical prior information
  for (i in 1:nrow(prior1)) {
    for (j in 1:ncol(prior1)) {
      if (prior1[i, j] != 0) {
        prior2[i, j] <- prior1[i, j]
        prior2[j, i] <- prior1[i, j]
      }
    }
  }
  prior2 <- ifelse(prior2 != 0, 1, 0)

  # network construction

  # network construction(start here)
  score1_tp <- array(0, dim = c(boot_time, 1))
  score1_fp <- array(0, dim = c(boot_time, 1))
  score1_fn <- array(0, dim = c(boot_time, 1))
  score1_tn <- array(0, dim = c(boot_time, 1))

  score2_tp <- array(0, dim = c(boot_time, 1))
  score2_fp <- array(0, dim = c(boot_time, 1))
  score2_fn <- array(0, dim = c(boot_time, 1))
  score2_tn <- array(0, dim = c(boot_time, 1))

  score3_tp <- array(0, dim = c(boot_time, 1))
  score3_fp <- array(0, dim = c(boot_time, 1))
  score3_fn <- array(0, dim = c(boot_time, 1))
  score3_tn <- array(0, dim = c(boot_time, 1))

  scorer_tp <- array(0, dim = c(boot_time, 1))
  scorer_fp <- array(0, dim = c(boot_time, 1))
  scorer_fn <- array(0, dim = c(boot_time, 1))
  scorer_tn <- array(0, dim = c(boot_time, 1))

  n_edges1 <- array(0, dim = c(boot_time, 1))
  n_edges2 <- array(0, dim = c(boot_time, 1))
  n_edges3 <- array(0, dim = c(boot_time, 1))
  n_edgesr <- array(0, dim = c(boot_time, 1))
  # loop start here
  for (r in 1:boot_time) {

    # Split the data into two parts
    valid_n <- sample(1:nrow(data), (1 / fold) * nrow(data), replace = FALSE) # default fold=5
    valid_data <- data[valid_n, ]
    train_data <- data[-valid_n, ]

    # get the parameters for regulization from training data(glasso-with prior)
    pc <- cov(train_data)

    # range of rho should be (0,1) but according to save the time,set to (0,0.1) as tested
    rho <- seq(0.01, 1, length = 100)
    bic <- matrix(NA, 100, 100)
    kappa <- rho
    rho_m <- NULL
    g_result <- NULL
    p_off_d <- NULL
    u <- matrix(1, nrow(prior2), ncol(prior2))
    for (i in 1:100) {
      for (j in 1:i) {
        rho_m <- rho[i] * u - kappa[j] * prior2
        g_result <- glasso::glasso(pc, rho_m)
        p_off_d <- sum(g_result$wi != 0 & col(pc) < row(pc))
        bic[i, j] <- -2 * (g_result$loglik) + p_off_d * log(nrow(train_data))
      }
    }
    pos <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)
    rho <- rho[pos[1]]
    kappa <- kappa[pos[2]]
    rho_m <- rho * u - kappa * prior2

    # Network construction with directional & prior information (trainning data)
    g_result <- glasso::glasso(pc, rho = rho_m)
    sigma_matrix <- g_result$wi
    niter <- g_result$niter
    if (niter == 10000) {
      stop("ERROR: Algorithm does not converge.")
    }
    pcor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor_matrix[i, j] <- -sigma_matrix[i, j] / sqrt(sigma_matrix[i, i] * sigma_matrix[j, j])
      }
    }

    # loglikelihood value of the training data
    # loglik_train <- g_result$loglik

    t_edges <- pcor_matrix

    # cut off = 0.1
    t_net_rhoadjusted2 <- as.data.frame(ifelse(abs(t_edges) >= 0.1 & row(t_edges) != col(t_edges), t_edges, 0))
    colnames(t_net_rhoadjusted2) <- colnames(data)
    rownames(t_net_rhoadjusted2) <- colnames(data)

    t_net_rhoadjusted_d <- matrix(0, nrow = nrow(t_net_rhoadjusted2), ncol = ncol(t_net_rhoadjusted2))
    for (i in 1:ncol(t_net_rhoadjusted_d)) {
      for (j in 1:nrow(t_net_rhoadjusted_d)) {
        if (prior1[i, j] != 0 & prior1[j, i] == 0) {
          t_net_rhoadjusted_d[i, j] <- t_net_rhoadjusted2[i, j]
        }
        if (prior1[i, j] == 0 & prior1[j, i] == 0) {
          t_net_rhoadjusted_d[i, j] <- t_net_rhoadjusted2[i, j]
          t_net_rhoadjusted_d[j, i] <- t_net_rhoadjusted2[j, i]
        }
        if (prior1[i, j] != 0 & prior1[j, i] != 0) {
          t_net_rhoadjusted_d[i, j] <- t_net_rhoadjusted2[i, j]
          t_net_rhoadjusted_d[j, i] <- t_net_rhoadjusted2[j, i]
        }
      }
    }
    colnames(t_net_rhoadjusted_d) <- colnames(data)
    rownames(t_net_rhoadjusted_d) <- colnames(data)

    # network with no prior information
    # get the regularization parameter
    rho <- seq(0.01, 1, length = 100)
    bic <- rho
    g_result <- NULL
    p_off_d <- NULL
    pc <- var(train_data)
    for (i in 1:100) {
      g_result <- glasso::glasso(pc, rho[i])
      p_off_d <- sum(g_result$wi != 0 & col(pc) < row(pc))
      bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(train_data))
    }
    best <- which.min(bic)
    rho <- rho[best]
    # get the network without prior
    g_result <- glasso::glasso(pc, rho = rho)
    sigma_matrix <- g_result$wi
    niter <- g_result$niter

    if (niter == 10000) {
      stop("ERROR: Algorithm does nor converge.")
    }

    pcor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor_matrix[i, j] <- -sigma_matrix[i, j] / sqrt(sigma_matrix[i, i] * sigma_matrix[j, j])
      }
    }

    t_edges_rho <- pcor_matrix
    t_edges_rho <- ifelse(abs(t_edges_rho) >= 0.1, t_edges_rho, 0)

    # network random
    eta_a <- sum(t_edges_rho != 0) / (ncol(train_data) * (ncol(train_data) - 1) * 2)
    random_pcor <- GeneNet::ggm.simulate.pcor(ncol(data), etaA = eta_a)

    # t_net in signaling way(0,1)
    t_net_rhoadjusted_d <- ifelse(col(t_net_rhoadjusted_d) != row(t_net_rhoadjusted_d) & t_net_rhoadjusted_d > 0, 1,
      ifelse(row(t_net_rhoadjusted_d) != col(t_net_rhoadjusted_d) & t_net_rhoadjusted_d < 0, -1, 0)
    )

    t_net_rhoadjusted2 <- ifelse(col(t_net_rhoadjusted2) != row(t_net_rhoadjusted2) & t_net_rhoadjusted2 > 0, 1,
      ifelse(row(t_net_rhoadjusted2) != col(t_net_rhoadjusted2) & t_net_rhoadjusted2 < 0, -1, 0)
    )

    t_net_rho <- ifelse(row(t_edges_rho) != col(t_edges_rho) & t_edges_rho > 0, 1,
      ifelse(row(t_edges_rho) != col(t_edges_rho) & t_edges_rho < 0, -1, 0)
    )

    t_net_r <- ifelse(col(random_pcor) != row(random_pcor) & random_pcor > 0, 1,
      ifelse(row(random_pcor) != col(random_pcor) & random_pcor < 0, -1, 0)
    )
    # number of non-zero edges predicted
    n_edges1[r] <- sum(t_net_rhoadjusted_d != 0)
    n_edges2[r] <- sum(t_net_rhoadjusted2 != 0)
    n_edges3[r] <- sum(t_net_rho != 0)
    n_edgesr[r] <- sum(t_net_r != 0)

    # number of zero edges predicted
    n0_edges1[r] <- sum(t_net_rhoadjusted_d == 0)
    n0_edges2[r] <- sum(t_net_rhoadjusted2 == 0)
    n0_edges3[r] <- sum(t_net_rho == 0)
    n0_edgesr[r] <- sum(t_net_r == 0)

    # validation data network generation
    # get the parameters for regulization from valid data (glasso-with prior)
    pc <- cov(valid_data)

    # range of rho should be (0,1) but according to save the time,set to (0,0.1) as tested
    rho <- seq(0.01, 1, length = 100)
    bic <- matrix(NA, 100, 100)
    kappa <- rho
    rho_m <- NULL
    g_result <- NULL
    p_off_d <- NULL
    u <- matrix(1, nrow(prior2), ncol(prior2))
    for (i in 1:100) {
      for (j in 1:i) {
        rho_m <- rho[i] * u - kappa[j] * prior2
        g_result <- glasso::glasso(pc, rho_m)
        p_off_d <- sum(g_result$wi != 0 & col(pc) < row(pc))
        bic[i, j] <- -2 * (g_result$loglik) + p_off_d * log(nrow(valid_data))
      }
    }
    pos <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)
    rho <- rho[pos[1]]
    kappa <- kappa[pos[2]]
    rho_m <- rho * u - kappa * prior2

    # Network construction with directional & prior information (valid_data)
    g_result <- glasso::glasso(pc, rho = rho_m)
    sigma_matrix <- g_result$wi
    niter <- g_result$niter

    if (niter == 10000) {
      stop("ERROR: Algorithm does nor converge.")
    }

    pcor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor_matrix[i, j] <- -sigma_matrix[i, j] / sqrt(sigma_matrix[i, i] * sigma_matrix[j, j])
      }
    }

    # loglikelihood value of valid data
    # loglik_valid <- g_result$loglik

    t_edges <- pcor_matrix

    # cut off = 0.1
    t_net_rhoadjusted2_valid <- as.data.frame(ifelse(abs(t_edges) >= 0.1 & row(t_edges) != col(t_edges), t_edges, 0))
    colnames(t_net_rhoadjusted2_valid) <- colnames(data)
    rownames(t_net_rhoadjusted2_valid) <- colnames(data)

    t_net_rhoadjusted_d_valid <- matrix(0, nrow = nrow(t_net_rhoadjusted2_valid), ncol = ncol(t_net_rhoadjusted2_valid))
    for (i in 1:ncol(t_net_rhoadjusted_d_valid)) {
      for (j in 1:nrow(t_net_rhoadjusted_d_valid)) {
        if (prior1[i, j] != 0 & prior1[j, i] == 0) {
          t_net_rhoadjusted_d_valid[i, j] <- t_net_rhoadjusted2_valid[i, j]
        }
        if (prior1[i, j] == 0 & prior1[j, i] == 0) {
          t_net_rhoadjusted_d_valid[i, j] <- t_net_rhoadjusted2_valid[i, j]
          t_net_rhoadjusted_d_valid[j, i] <- t_net_rhoadjusted2_valid[j, i]
        }
        if (prior1[i, j] != 0 & prior1[j, i] != 0) {
          t_net_rhoadjusted_d_valid[i, j] <- t_net_rhoadjusted2_valid[i, j]
          t_net_rhoadjusted_d_valid[j, i] <- t_net_rhoadjusted2_valid[j, i]
        }
      }
    }
    colnames(t_net_rhoadjusted_d_valid) <- colnames(data)
    rownames(t_net_rhoadjusted_d_valid) <- colnames(data)

    # network without prior without direction glasso valid data
    rho <- seq(0.01, 1, length = 100)
    bic <- rho
    g_result <- NULL
    p_off_d <- NULL
    pc <- var(valid_data)
    for (i in 1:100) {
      g_result <- glasso::glasso(pc, rho[i])
      p_off_d <- sum(g_result$wi != 0 & col(pc) < row(pc))
      bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(valid_data))
    }
    best <- which.min(bic)
    rho <- rho[best]
    # get the network without prior
    g_result <- glasso::glasso(pc, rho = rho)
    sigma_matrix <- g_result$wi
    pcor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor_matrix[i, j] <- -sigma_matrix[i, j] / sqrt(sigma_matrix[i, i] * sigma_matrix[j, j])
      }
    }

    t_edges_rho_valid <- pcor_matrix
    t_edges_rho_valid <- ifelse(abs(t_edges_rho_valid) >= 0.1, t_edges_rho_valid, 0)

    # random_network(valid_data)

    random_pcor_valid <- GeneNet::ggm.simulate.pcor(ncol(data), etaA = sum(t_edges_rho_valid != 0) /
      (ncol(valid_data) * (ncol(valid_data) - 1) * 2))

    # t_net in signaling way(0,1)
    t_net_rhoadjusted_d_valid <- ifelse(col(t_net_rhoadjusted_d_valid) != row(t_net_rhoadjusted_d_valid) &
      t_net_rhoadjusted_d_valid > 0, 1,
    ifelse(row(t_net_rhoadjusted_d_valid) != col(t_net_rhoadjusted_d_valid) & t_net_rhoadjusted_d_valid < 0, -1, 0)
    )

    t_net_rhoadjusted2_valid <- ifelse(col(t_net_rhoadjusted2_valid) != row(t_net_rhoadjusted2_valid) &
      t_net_rhoadjusted2_valid > 0, 1,
    ifelse(row(t_net_rhoadjusted2_valid) != col(t_net_rhoadjusted2_valid) & t_net_rhoadjusted2_valid < 0, -1, 0)
    )

    t_net_rho_valid <- ifelse(row(t_edges_rho_valid) != col(t_edges_rho_valid) & t_edges_rho_valid > 0, 1,
      ifelse(row(t_edges_rho_valid) != col(t_edges_rho_valid) & t_edges_rho_valid < 0, -1, 0)
    )

    t_net_r_valid <- ifelse(col(random_pcor_valid) != row(random_pcor_valid) & random_pcor_valid > 0, 1,
      ifelse(row(random_pcor_valid) != col(random_pcor_valid) & random_pcor_valid < 0, -1, 0)
    )


    # 5 fold cross validation as comparisn of the trainning network with the signal of the valid data
    score1_1 <- 0
    score1_2 <- 0
    score1_3 <- 0
    score1_4 <- 0

    score2_1 <- 0
    score2_2 <- 0
    score2_3 <- 0
    score2_4 <- 0

    score3_1 <- 0
    score3_2 <- 0
    score3_3 <- 0
    score3_4 <- 0

    score_r_1 <- 0
    score_r_2 <- 0
    score_r_3 <- 0
    score_r_4 <- 0

    # True Positive#
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t_net_rhoadjusted_d_valid[t, p] == t_net_rhoadjusted_d[t, p] & t_net_rhoadjusted_d_valid[t, p] != 0) {
          score1_1 <- score1_1 + 1
        }
        if (t_net_rhoadjusted2_valid[t, p] == t_net_rhoadjusted2[t, p] & t_net_rhoadjusted2_valid[t, p] != 0) {
          score2_1 <- score2_1 + 1
        }
        if (t_net_rho_valid[t, p] == t_net_rho[t, p] & t_net_rho_valid[t, p] != 0) {
          score3_1 <- score3_1 + 1
        }
        if (t_net_r_valid[t, p] == t_net_r[t, p] & t_net_r_valid[t, p] != 0) {
          score_r_1 <- score_r_1 + 1
        }
      }
    }

    # False Positive
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t_net_rhoadjusted_d_valid[t, p] != t_net_rhoadjusted_d[t, p] &
          t_net_rhoadjusted_d[t, p] != 0 & t_net_rhoadjusted_d_valid[t, p] != 0) {
          score1_2 <- score1_2 + 1
        }
        if (t_net_rhoadjusted2_valid[t, p] != t_net_rhoadjusted2[t, p] &
          t_net_rhoadjusted2[t, p] != 0 & t_net_rhoadjusted2_valid[t, p] != 0) {
          score2_2 <- score2_2 + 1
        }
        if (t_net_rho_valid[t, p] != t_net_rho[t, p] & t_net_rho[t, p] != 0 &
          t_net_rho_valid[t, p] != 0) {
          score3_2 <- score3_2 + 1
        }
        if (t_net_r_valid[t, p] != t_net_r[t, p] & t_net_r[t, p] != 0 &
          t_net_r_valid[t, p] != 0) {
          score_r_2 <- score_r_2 + 1
        }
      }
    }

    # False Negative
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t_net_rhoadjusted_d_valid[t, p] != t_net_rhoadjusted_d[t, p] & t_net_rhoadjusted_d[t, p] == 0) {
          score1_3 <- score1_3 + 1
        }
        if (t_net_rhoadjusted2_valid[t, p] != t_net_rhoadjusted2[t, p] & t_net_rhoadjusted2[t, p] == 0) {
          score2_3 <- score2_3 + 1
        }
        if (t_net_rho_valid[t, p] != t_net_rho[t, p] & t_net_rho[t, p] == 0) {
          score3_3 <- score3_3 + 1
        }
        if (t_net_r_valid[t, p] != t_net_r[t, p] & t_net_r[t, p] == 0) {
          score_r_3 <- score_r_3 + 1
        }
      }
    }

    # True Negative
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t_net_rhoadjusted_d_valid[t, p] == t_net_rhoadjusted_d[t, p] & t_net_rhoadjusted_d_valid[t, p] == 0) {
          score1_4 <- score1_4 + 1
        }
        if (t_net_rhoadjusted2_valid[t, p] == t_net_rhoadjusted2[t, p] & t_net_rhoadjusted2_valid[t, p] == 0) {
          score2_4 <- score2_4 + 1
        }
        if (t_net_rho_valid[t, p] == t_net_rho[t, p] & t_net_rho_valid[t, p] == 0) {
          score3_4 <- score3_4 + 1
        }
        if (t_net_r_valid[t, p] == t_net_r[t, p] & t_net_r_valid[t, p] == 0) {
          score_r_4 <- score_r_4 + 1
        }
      }
    }
    # score in table###and the explaination##3

    score1_tp[r] <- score1_1
    score1_fp[r] <- score1_2
    score1_fn[r] <- score1_3
    score1_tn[r] <- score1_4

    score2_tp[r] <- score2_1
    score2_fp[r] <- score2_2
    score2_fn[r] <- score2_3
    score2_tn[r] <- score2_4

    score3_tp[r] <- score3_1
    score3_fp[r] <- score3_2
    score3_fn[r] <- score3_3
    score3_tn[r] <- score3_4

    scorer_tp[r] <- score_r_1
    scorer_fp[r] <- score_r_2
    scorer_fn[r] <- score_r_3
    scorer_tn[r] <- score_r_4
  }
  # loop end here,write out the result here with no overlap
  # true positive
  score1_p <- mean(score1_tp / n_edges1)
  score_b_1 <- score1_fn / n_edges1
  # false negative
  score1_f <- mean(score1_fn / n0_edges1)
  score0_b_1 <- score1_fn / n0_edges1

  # true positive
  score2_p <- mean(score2_tp / n_edges2)
  score_b_2 <- score2_tp / n_edges2
  # false negative
  score2_f <- mean(score2_fn / n0_edges2)
  score0_b_2 <- score1_fn / n0_edges2

  # true positive
  score3_p <- mean(score3_tp / n_edges3)
  score_b_3 <- score3_tp / n_edges3
  # false negative
  score3_f <- mean(score1_fn / n0_edges3)
  score0_b_3 <- score1_fn / n0_edges3

  # true positive
  scorer_p <- mean(scorer_tp / n_edgesr)
  score_b_r <- scorer_tp / n_edgesr
  # false negative
  scorer_f <- mean(scorer_fn / n0_edgesr)
  score0_b_r <- scorer_fn / n0_edgesr

  result <- list(
    score1_p = score1_p, score_b_1 = score_b_1, score1_f = score1_f, score0_b_1 = score0_b_1,
    score2_p = score2_p, score_b_2 = score_b_2, score2_f = score2_f, score0_b_2 = score0_b_2,
    score3_p = score3_p, score_b_3 = score_b_3, score3_f = score3_f, score0_b_3 = score0_b_3,
    scorer_p = scorer_p, score_b_r = score_b_r, scorer_f = scorer_f, score0_b_r = score0_b_r
  )

  return(result)
}
