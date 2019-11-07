#' Run k-fold cross validation for data.(reference from CVglasso)
#'
#'  @param data option to provide a nxp data matrix. Each row corresponds to a single observation
#'   and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is
#'  \code{NULL} and \code{data} is provided instead then \code{S} will be computed automatically.
#' @param rho positive tuning parameter for elastic net penalty. Default at seq(0.01,1,length=100)
#' @param kappa positive scaler parameter for biology prior contribution. Default at seq(0.01,1,length=100)
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix
#' (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param tol convergence tolerance. Iterations will stop when the average absolute difference
#'  in parameter estimates in less than \code{tol} times multiple. Defaults to 1e-4.
#' @param maxit maximum number of iterations. Defaults to 1e4.
#' @param adjmaxit adjusted maximum number of iterations. During cross validation this option
#' allows the user to adjust the maximum number of iterations after the first \code{lam} tuning
#'  parameter has converged. This option is intended to be paired with \code{warm} starts and
#'  allows for 'one-step' estimators. Defaults to NULL.
#' @param k_fold specify the number of folds for cross validation.
#' @param crit_cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param algorithm
#'
#' @return returns list of returns which includes:
#' \item{lam}{optimized penalty parameter through traing data.}
#' \item{avg_error}{average cross validation error across all folds.}
#' \item{cv_error}{cross validation errors.}
#'
#' @importFrom glasso glasso
#'
#' @concept zeptosensPkg
#' @export

cv_glasso <- function(data = NULL,
                      s_matrix = NULL,
                      prior,
                      rho = seq(0.01, 1, length = 100),
                      kappa = seq(0.01, 1, length = 100),
                      diagonal = FALSE,
                      tol = 1e-04,
                      maxit = 10000,
                      adjmaxit = NULL,
                      k_fold = 2,
                      crit_cv = c("loglik", "AIC", "BIC"),
                      start = c("warm", "cold"),
                      algorithm = c("data_driven", "hybrid_driven")) {

  # match values
  crit_cv <- match.arg(crit_cv)
  algorithm <- match.arg(algorithm)
  start <- match.arg(start)

  # initialize
  initmaxit <- maxit
  s_train <- s_matrix
  s_valid <- s_matrix
  crossvalid_errors <- array(0, k_fold)

  # no need to create folds if k_fold <- 1
  if (k_fold == 1) {

    # set sample size
    n <- nrow(s_matrix)
  } else {

    # designate folds and shuffle -- ensures randomized folds
    n <- nrow(data)
    ind <- sample(n)
  }

  # parse data into folds and perform CV
  for (k in 1:k_fold) {
    if (k_fold > 1) {

      # training set
      leave_out <- ind[(1 + floor((k - 1) * n / k_fold)):floor(k *
        n / k_fold)]
      data_train <- data[-leave_out, , drop = FALSE]
      data_bar <- apply(data_train, 2, mean)
      data_train <- scale(data_train, center = data_bar, scale = FALSE)

      # validation set
      data_valid <- data[leave_out, , drop = FALSE]
      data_valid <- scale(data_valid, center = data_bar, scale = FALSE)
      n <- nrow(data_valid)

      # sample covariances
      s_train <- crossprod(data_train) / (dim(data_train)[1])
      s_valid <- crossprod(data_valid) / (dim(data_valid)[1])
    }

    # re-initialize values for each fold
    maxit <- initmaxit
    init <- s_train
    init_omega <- diag(ncol(s_train))

    # for tuning parameters
    # choosing optimal parameter tuning parameter for training data
    if (algorithm == "data_driven") {
      optimize_param <- zeptosensPkg::optimize_parameter_dat(
        data = data_train,
        rho = rho
      )
      lam_ <- optimize_param$rho
    }

    if (algorithm == "hybrid_driven") {
      optimize_param <- zeptosensPkg::optimize_parameter_hybrid(
        data = data_train,
        rho = rho,
        kappa = kappa,
        prior = prior
      )
      lam_ <- optimize_param$rho_m
    }

    # initial sigma
    if (!diagonal) {

      # provide estimate that is pd and dual feasible
      s_minus <- s_train
      diag(s_minus) <- 0
      alpha <- min(c(min(rho) / max(abs(s_minus)), 1))
      init <- (1 - alpha) * s_train
      diag(init) <- diag(s_train)
    }

    # update diagonal elements of init, if necessary
    if (diagonal) {
      diag(init) <- diag(s_train) + lam_
    }

    # compute the penalized likelihood precision matrix
    # estimator
    glasso_result <- glasso::glasso(
      s = s_train, rho = lam_, thr = tol,
      maxit = maxit, penalize.diagonal = diagonal,
      start = "warm", w.init = init, wi.init = init_omega
    )

    if (start == "warm") {

      # option to save initial values for warm starts
      init <- glasso_result$w
      init_omega <- glasso_result$wi
      maxit <- adjmaxit
    }

    # compute the observed negative validation loglikelihood
    # (close enoug)
    crossvalid_errors[k] <- (nrow(data) / 2) * (sum(glasso_result$wi *
      s_valid) - determinant(glasso_result$wi, logarithm = TRUE)$modulus[1])

    # update for crit_cv, if necessary
    if (crit_cv == "AIC") {
      crossvalid_errors[k] <- crossvalid_errors[k] + sum(glasso_result$wi != 0)
    }
    if (crit_cv == "BIC") {
      crossvalid_errors[k] <- crossvalid_errors[k] + sum(glasso_result$wi != 0) * log(nrow(data)) / 2
    }
  }

  # determine optimal tuning parameters
  avg_error <- mean(crossvalid_errors)


  # return best lam and alpha values
  return(list(lam = lam_, avg_error = avg_error, cv_error = crossvalid_errors))
}
