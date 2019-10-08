#' Run Cross validation for data through the chosen algorithm.
#'
#' @param data input expression data. Coloumns as the gene, rows as the sample.With colnames as the gene tags, rownames as the sample tags.
#' @param prior prior information matrix with colnames and rownames as the gene tags.
#' @param boot.time Bootstrap time mannually set.Default at 1000.
#' @param fold The fold for training and test dataset. Default at 5.
#' @return result of validation score. for random network and the network predicted with the algorithm.
#' @concept zeptosensPkg
#' @export
runCrossValidation <- function(data, prior, boot.time = 1000, cut.off = 0.1, fold = 5) {
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
  score1.tp <- array(0, dim = c(boot.time, 1))
  score1.fp <- array(0, dim = c(boot.time, 1))
  score1.fn <- array(0, dim = c(boot.time, 1))
  score1.tn <- array(0, dim = c(boot.time, 1))

  score2.tp <- array(0, dim = c(boot.time, 1))
  score2.fp <- array(0, dim = c(boot.time, 1))
  score2.fn <- array(0, dim = c(boot.time, 1))
  score2.tn <- array(0, dim = c(boot.time, 1))

  score3.tp <- array(0, dim = c(boot.time, 1))
  score3.fp <- array(0, dim = c(boot.time, 1))
  score3.fn <- array(0, dim = c(boot.time, 1))
  score3.tn <- array(0, dim = c(boot.time, 1))

  scorer.tp <- array(0, dim = c(boot.time, 1))
  scorer.fp <- array(0, dim = c(boot.time, 1))
  scorer.fn <- array(0, dim = c(boot.time, 1))
  scorer.tn <- array(0, dim = c(boot.time, 1))

  n_edges1 <- array(0, dim = c(boot.time, 1))
  n_edges2 <- array(0, dim = c(boot.time, 1))
  n_edges3 <- array(0, dim = c(boot.time, 1))
  n_edgesr <- array(0, dim = c(boot.time, 1))
  # loop start here
  for (r in 1:boot.time) {

    # Split the data into two parts
    valid_n <- sample(1:nrow(data), (1 / fold) * nrow(data), replace = F) # default fold=5
    valid_data <- data[valid_n, ]
    train_data <- data[-valid_n, ]

    # get the parameters for regulization from training data(glasso-with prior)
    pc <- cov(train_data)
    rho <- seq(0.01, 1, length = 100) # range of rho should be (0,1) but according to save the time,set to (0,0.1) as tested
    bic <- matrix(NA, 100, 100)
    kappa <- rho
    rho_m <- c()
    g.result <- c()
    p_off_d <- c()
    U <- matrix(1, nrow(prior2), ncol(prior2))
    for (i in 1:100) {
      for (j in 1:i) {
        rho_m <- rho[i] * U - kappa[j] * prior2
        g.result <- glasso(pc, rho_m)
        p_off_d <- sum(g.result$wi != 0 & col(pc) < row(pc))
        bic[i, j] <- -2 * (g.result$loglik) + p_off_d * log(nrow(train_data))
      }
    }
    pos <- which(bic == min(bic, na.rm = TRUE), arr.ind = T)
    rho <- rho[pos[1]]
    kappa <- kappa[pos[2]]
    rho_m <- rho * U - kappa * prior2

    # Network construction with directional & prior information (trainning data)
    g_result <- glasso(pc, rho = rho_m)
    sigma.matrix <- g_result$wi
    niter <- g_result$niter
    if (niter == 10000) {
      stop("ERROR: Algorithm does not converge.")
    }
    pcor.matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor.matrix[i, j] <- -sigma.matrix[i, j] / sqrt(sigma.matrix[i, i] * sigma.matrix[j, j])
      }
    }

    # loglikelihood value of the training data
    loglik_train <- g_result$loglik

    t.edges <- pcor.matrix

    # cut off = 0.1
    t.net.rhoadjusted2 <- as.data.frame(ifelse(abs(t.edges) >= 0.1 & row(t.edges) != col(t.edges), t.edges, 0))
    colnames(t.net.rhoadjusted2) <- colnames(data)
    rownames(t.net.rhoadjusted2) <- colnames(data)

    t.net.rhoadjusted.d <- matrix(0, nrow = nrow(t.net.rhoadjusted2), ncol = ncol(t.net.rhoadjusted2))
    for (i in 1:ncol(t.net.rhoadjusted.d)) {
      for (j in 1:nrow(t.net.rhoadjusted.d)) {
        if (prior1[i, j] != 0 & prior1[j, i] == 0) {
          t.net.rhoadjusted.d[i, j] <- t.net.rhoadjusted2[i, j]
        }
        if (prior1[i, j] == 0 & prior1[j, i] == 0) {
          t.net.rhoadjusted.d[i, j] <- t.net.rhoadjusted2[i, j]
          t.net.rhoadjusted.d[j, i] <- t.net.rhoadjusted2[j, i]
        }
        if (prior1[i, j] != 0 & prior1[j, i] != 0) {
          t.net.rhoadjusted.d[i, j] <- t.net.rhoadjusted2[i, j]
          t.net.rhoadjusted.d[j, i] <- t.net.rhoadjusted2[j, i]
        }
      }
    }
    colnames(t.net.rhoadjusted.d) <- colnames(data)
    rownames(t.net.rhoadjusted.d) <- colnames(data)

    # network with no prior information
    # get the regularization parameter
    rho <- seq(0.01, 1, length = 100)
    bic <- rho
    g_result <- c()
    p_off_d <- c()
    pc <- var(train_data)
    for (i in 1:100) {
      g_result <- glasso(pc, rho[i])
      p_off_d <- sum(g_result$wi != 0 & col(pc) < row(pc))
      bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(train_data))
    }
    best <- which.min(bic)
    rho <- rho[best]
    # get the network without prior
    g_result <- glasso(pc, rho = rho)
    sigma.matrix <- g_result$wi
    niter <- g_result$niter

    if (niter == 10000) {
      stop("ERROR: Algorithm does nor converge.")
    }

    pcor.matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor.matrix[i, j] <- -sigma.matrix[i, j] / sqrt(sigma.matrix[i, i] * sigma.matrix[j, j])
      }
    }

    t.edges.rho <- pcor.matrix
    t.edges.rho <- ifelse(abs(t.edges.rho) >= 0.1, t.edges.rho, 0)

    # network random
    random.pcor <- ggm.simulate.pcor(ncol(data), etaA = sum(t.edges.rho != 0) / (ncol(train_data) * (ncol(train_data) - 1) * 2))

    # t.net in signaling way(0,1)
    t.net.rhoadjusted.d <- ifelse(col(t.net.rhoadjusted.d) != row(t.net.rhoadjusted.d) & t.net.rhoadjusted.d > 0, 1,
      ifelse(row(t.net.rhoadjusted.d) != col(t.net.rhoadjusted.d) & t.net.rhoadjusted.d < 0, -1, 0)
    )

    t.net.rhoadjusted2 <- ifelse(col(t.net.rhoadjusted2) != row(t.net.rhoadjusted2) & t.net.rhoadjusted2 > 0, 1,
      ifelse(row(t.net.rhoadjusted2) != col(t.net.rhoadjusted2) & t.net.rhoadjusted2 < 0, -1, 0)
    )

    t.net.rho <- ifelse(row(t.edges.rho) != col(t.edges.rho) & t.edges.rho > 0, 1,
      ifelse(row(t.edges.rho) != col(t.edges.rho) & t.edges.rho < 0, -1, 0)
    )

    t.net.r <- ifelse(col(random.pcor) != row(random.pcor) & random.pcor > 0, 1,
      ifelse(row(random.pcor) != col(random.pcor) & random.pcor < 0, -1, 0)
    )
    # number of non-zero edges predicted
    n_edges1[r] <- sum(t.net.rhoadjusted.d != 0)
    n_edges2[r] <- sum(t.net.rhoadjusted2 != 0)
    n_edges3[r] <- sum(t.net.rho != 0)
    n_edgesr[r] <- sum(t.net.r != 0)

    # number of zero edges predicted
    n0_edges1[r] <- sum(t.net.rhoadjusted.d == 0)
    n0_edges2[r] <- sum(t.net.rhoadjusted2 == 0)
    n0_edges3[r] <- sum(t.net.rho == 0)
    n0_edgesr[r] <- sum(t.net.r == 0)

    # validation data network generation
    # get the parameters for regulization from valid data (glasso-with prior)
    pc <- cov(valid_data)
    rho <- seq(0.01, 1, length = 100) # range of rho should be (0,1) but according to save the time,set to (0,0.1) as tested
    bic <- matrix(NA, 100, 100)
    kappa <- rho
    rho_m <- c()
    g.result <- c()
    p_off_d <- c()
    U <- matrix(1, nrow(prior2), ncol(prior2))
    for (i in 1:100) {
      for (j in 1:i) {
        rho_m <- rho[i] * U - kappa[j] * prior2
        g.result <- glasso(pc, rho_m)
        p_off_d <- sum(g.result$wi != 0 & col(pc) < row(pc))
        bic[i, j] <- -2 * (g.result$loglik) + p_off_d * log(nrow(valid_data))
      }
    }
    pos <- which(bic == min(bic, na.rm = TRUE), arr.ind = T)
    rho <- rho[pos[1]]
    kappa <- kappa[pos[2]]
    rho_m <- rho * U - kappa * prior2

    # Network construction with directional & prior information (valid_data)
    g_result <- glasso(pc, rho = rho_m)
    sigma.matrix <- g_result$wi
    niter <- g_result$niter

    if (niter == 10000) {
      stop("ERROR: Algorithm does nor converge.")
    }

    pcor.matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor.matrix[i, j] <- -sigma.matrix[i, j] / sqrt(sigma.matrix[i, i] * sigma.matrix[j, j])
      }
    }

    # loglikelihood value of valid data
    loglik_valid <- g_result$loglik


    t.edges <- pcor.matrix

    # cut off = 0.1
    t.net.rhoadjusted2_valid <- as.data.frame(ifelse(abs(t.edges) >= 0.1 & row(t.edges) != col(t.edges), t.edges, 0))
    colnames(t.net.rhoadjusted2_valid) <- colnames(data)
    rownames(t.net.rhoadjusted2_valid) <- colnames(data)

    t.net.rhoadjusted.d_valid <- matrix(0, nrow = nrow(t.net.rhoadjusted2_valid), ncol = ncol(t.net.rhoadjusted2_valid))
    for (i in 1:ncol(t.net.rhoadjusted.d_valid)) {
      for (j in 1:nrow(t.net.rhoadjusted.d_valid)) {
        if (prior1[i, j] != 0 & prior1[j, i] == 0) {
          t.net.rhoadjusted.d_valid[i, j] <- t.net.rhoadjusted2_valid[i, j]
        }
        if (prior1[i, j] == 0 & prior1[j, i] == 0) {
          t.net.rhoadjusted.d_valid[i, j] <- t.net.rhoadjusted2_valid[i, j]
          t.net.rhoadjusted.d_valid[j, i] <- t.net.rhoadjusted2_valid[j, i]
        }
        if (prior1[i, j] != 0 & prior1[j, i] != 0) {
          t.net.rhoadjusted.d_valid[i, j] <- t.net.rhoadjusted2_valid[i, j]
          t.net.rhoadjusted.d_valid[j, i] <- t.net.rhoadjusted2_valid[j, i]
        }
      }
    }
    colnames(t.net.rhoadjusted.d_valid) <- colnames(data)
    rownames(t.net.rhoadjusted.d_valid) <- colnames(data)

    # network without prior without direction glasso valid data
    rho <- seq(0.01, 1, length = 100)
    bic <- rho
    g_result <- c()
    p_off_d <- c()
    pc <- var(valid_data)
    for (i in 1:100) {
      g_result <- glasso(pc, rho[i])
      p_off_d <- sum(g_result$wi != 0 & col(pc) < row(pc))
      bic[i] <- -2 * (g_result$loglik) + p_off_d * log(nrow(valid_data))
    }
    best <- which.min(bic)
    rho <- rho[best]
    # get the network without prior
    g_result <- glasso(pc, rho = rho)
    sigma.matrix <- g_result$wi
    pcor.matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      for (j in 1:ncol(data)) {
        pcor.matrix[i, j] <- -sigma.matrix[i, j] / sqrt(sigma.matrix[i, i] * sigma.matrix[j, j])
      }
    }

    t.edges.rho_valid <- pcor.matrix
    t.edges.rho_valid <- ifelse(abs(t.edges.rho_valid) >= 0.1, t.edges.rho_valid, 0)

    # random_network(valid_data)

    random.pcor_valid <- ggm.simulate.pcor(ncol(data), etaA = sum(t.edges.rho_valid != 0) / (ncol(valid_data) * (ncol(valid_data) - 1) * 2))

    # t.net in signaling way(0,1)
    t.net.rhoadjusted.d_valid <- ifelse(col(t.net.rhoadjusted.d_valid) != row(t.net.rhoadjusted.d_valid) & t.net.rhoadjusted.d_valid > 0, 1,
      ifelse(row(t.net.rhoadjusted.d_valid) != col(t.net.rhoadjusted.d_valid) & t.net.rhoadjusted.d_valid < 0, -1, 0)
    )

    t.net.rhoadjusted2_valid <- ifelse(col(t.net.rhoadjusted2_valid) != row(t.net.rhoadjusted2_valid) & t.net.rhoadjusted2_valid > 0, 1,
      ifelse(row(t.net.rhoadjusted2_valid) != col(t.net.rhoadjusted2_valid) & t.net.rhoadjusted2_valid < 0, -1, 0)
    )

    t.net.rho_valid <- ifelse(row(t.edges.rho_valid) != col(t.edges.rho_valid) & t.edges.rho_valid > 0, 1,
      ifelse(row(t.edges.rho_valid) != col(t.edges.rho_valid) & t.edges.rho_valid < 0, -1, 0)
    )

    t.net.r_valid <- ifelse(col(random.pcor_valid) != row(random.pcor_valid) & random.pcor_valid > 0, 1,
      ifelse(row(random.pcor_valid) != col(random.pcor_valid) & random.pcor_valid < 0, -1, 0)
    )


    # 5 fold cross validation as comparisn of the trainning network with the signal of the valid data
    score1.1 <- 0
    score1.2 <- 0
    score1.3 <- 0
    score1.4 <- 0

    score2.1 <- 0
    score2.2 <- 0
    score2.3 <- 0
    score2.4 <- 0

    score3.1 <- 0
    score3.2 <- 0
    score3.3 <- 0
    score3.4 <- 0

    score_r.1 <- 0
    score_r.2 <- 0
    score_r.3 <- 0
    score_r.4 <- 0

    # True Positive#
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t.net.rhoadjusted.d_valid[t, p] == t.net.rhoadjusted.d[t, p] & t.net.rhoadjusted.d_valid[t, p] != 0) {
          score1.1 <- score1.1 + 1
        }
        if (t.net.rhoadjusted2_valid[t, p] == t.net.rhoadjusted2[t, p] & t.net.rhoadjusted2_valid[t, p] != 0) {
          score2.1 <- score2.1 + 1
        }
        if (t.net.rho_valid[t, p] == t.net.rho[t, p] & t.net.rho_valid[t, p] != 0) {
          score3.1 <- score3.1 + 1
        }
        if (t.net.r_valid[t, p] == t.net.r[t, p] & t.net.r_valid[t, p] != 0) {
          score_r.1 <- score_r.1 + 1
        }
      }
    }

    # False Positive
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t.net.rhoadjusted.d_valid[t, p] != t.net.rhoadjusted.d[t, p] & t.net.rhoadjusted.d[t, p] != 0 & t.net.rhoadjusted.d_valid[t, p] != 0) {
          score1.2 <- score1.2 + 1
        }
        if (t.net.rhoadjusted2_valid[t, p] != t.net.rhoadjusted2[t, p] & t.net.rhoadjusted2[t, p] != 0 & t.net.rhoadjusted2_valid[t, p] != 0) {
          score2.2 <- score2.2 + 1
        }
        if (t.net.rho_valid[t, p] != t.net.rho[t, p] & t.net.rho[t, p] != 0 & t.net.rho_valid[t, p] != 0) {
          score3.2 <- score3.2 + 1
        }
        if (t.net.r_valid[t, p] != t.net.r[t, p] & t.net.r[t, p] != 0 & t.net.r_valid[t, p] != 0) {
          score_r.2 <- score_r.2 + 1
        }
      }
    }

    # False Negative#
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t.net.rhoadjusted.d_valid[t, p] != t.net.rhoadjusted.d[t, p] & t.net.rhoadjusted.d[t, p] == 0) {
          score1.3 <- score1.3 + 1
        }
        if (t.net.rhoadjusted2_valid[t, p] != t.net.rhoadjusted2[t, p] & t.net.rhoadjusted2[t, p] == 0) {
          score2.3 <- score2.3 + 1
        }
        if (t.net.rho_valid[t, p] != t.net.rho[t, p] & t.net.rho[t, p] == 0) {
          score3.3 <- score3.3 + 1
        }
        if (t.net.r_valid[t, p] != t.net.r[t, p] & t.net.r[t, p] == 0) {
          score_r.3 <- score_r.3 + 1
        }
      }
    }

    # True Negative#
    for (t in 1:ncol(data)) {
      for (p in 1:ncol(data)) {
        if (t.net.rhoadjusted.d_valid[t, p] == t.net.rhoadjusted.d[t, p] & t.net.rhoadjusted.d_valid[t, p] == 0) {
          score1.4 <- score1.4 + 1
        }
        if (t.net.rhoadjusted2_valid[t, p] == t.net.rhoadjusted2[t, p] & t.net.rhoadjusted2_valid[t, p] == 0) {
          score2.4 <- score2.4 + 1
        }
        if (t.net.rho_valid[t, p] == t.net.rho[t, p] & t.net.rho_valid[t, p] == 0) {
          score3.4 <- score3.4 + 1
        }
        if (t.net.r_valid[t, p] == t.net.r[t, p] & t.net.r_valid[t, p] == 0) {
          score_r.4 <- score_r.4 + 1
        }
      }
    }
    # score in table###and the explaination##3

    score1.tp[r] <- score1.1
    score1.fp[r] <- score1.2
    score1.fn[r] <- score1.3
    score1.tn[r] <- score1.4

    score2.tp[r] <- score2.1
    score2.fp[r] <- score2.2
    score2.fn[r] <- score2.3
    score2.tn[r] <- score2.4

    score3.tp[r] <- score3.1
    score3.fp[r] <- score3.2
    score3.fn[r] <- score3.3
    score3.tn[r] <- score3.4

    scorer.tp[r] <- score_r.1
    scorer.fp[r] <- score_r.2
    scorer.fn[r] <- score_r.3
    scorer.tn[r] <- score_r.4
  }
  # loop end here,write out the result here with no overlap
  # true positive
  score1.p <- mean(score1.tp / n_edges1)
  score.b.1 <- score1.fn / n_edges1
  # false negative
  score1.f <- mean(score1.fn / n0_edges1)
  score0.b.1 <- score1.fn / n0_edges1

  # true positive
  score2.p <- mean(score2.tp / n_edges2)
  score.b.2 <- score2.tp / n_edges2
  # false negative
  score2.f <- mean(score2.fn / n0_edges2)
  score0.b.2 <- score1.fn / n0_edges2

  # true positive
  score3.p <- mean(score3.tp / n_edges3)
  score.b.3 <- score3.tp / n_edges3
  # false negative
  score3.f <- mean(score1.fn / n0_edges3)
  score0.b.3 <- score1.fn / n0_edges3

  # true positive
  scorer.p <- mean(scorer.tp / n_edgesr)
  score.b.r <- scorer.tp / n_edgesr
  # false negative
  scorer.f <- mean(scorer.fn / n0_edgesr)
  score0.b.r <- scorer.fn / n0_edgesr

  result <- list(
    score1.p = score1.p, score.b.1 = score.b.1, score1.f = score1.f, score0.b.1 = score0.b.1,
    score2.p = score2.p, score.b.2 = score.b.2, score2.f = score2.f, score0.b.2 = score0.b.2,
    score3.p = score3.p, score.b.3 = score.b.3, score3.f = score3.f, score0.b.3 = score0.b.3,
    scorer.p = scorer.p, score.b.r = score.b.r, scorer.f = scorer.f, score0.b.r = score0.b.r
  )
  return(result)
}
