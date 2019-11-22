#' Run cross validation through data.
#'
#' @param n_prot antibody number of input data.
#' @param proteomic_responses input drug perturbation data. With columns as antibody, rows as samples.
#' @param max_dist maximum distance between two antibody. (Default at 1)
#' @param mab_to_genes A list of antibodies, their associated genes, modification sites and effect.
#' @param dist_file a distance file an edgelist with a third column which is the network distance
#' between the genes in the interaction
#' @param verbose whether to show debugging information
#' @param data  input proteomics dataset for network inference. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param prior prior information matrix of gene interaction, with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags.Can be inferred from Public data source (for example:SignedPC).
#' @param algorithm flexible toolbox implementing network estimating algorithms for robustness test.
#' (\code{data_driven},or \code{hybrid_driven}).
#' @param boot_time bootstrap time mannually set.Default at 100.
#' @param crit_cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param k_fold specify the number of folds for cross validation.(Default at 5)
#'
#' @return returns list of cross validation results which includes:
#' \item{cv_error}{cross validation errors of data.}
#' \item{cv_error_r}{cross validation errors of random data.}
#' \item{cv_error_diff}{cross validation errors difference between random data and data.}
#'
#' @importFrom glasso glasso
#'
#' @concept zeptosensPkg
#' @export

run_cross_validation <- function(data,
                                 prior = NULL,
                                 n_prot,
                                 k_fold = 5,
                                 proteomic_responses,
                                 mab_to_genes,
                                 dist_file = NULL,
                                 max_dist = 1,
                                 verbose = FALSE,
                                 boot_time = 100,
                                 crit_cv = "loglik",
                                 algorithm = c("data_driven", "hybrid_driven")) {
  # Match the data
  index <- colnames(data[which(colnames(data) %in% colnames(proteomic_responses))])
  data[, index]
  proteomic_responses <- proteomic_responses[, index]

  # Import prior biology information/NULL from SignedPC
  if (is.null(prior)) {
    network_ref <- zeptosensPkg::predict_bio_network(
      n_prot = n_prot,
      proteomic_responses = proteomic_responses,
      max_dist = max_dist,
      mab_to_genes = mab_to_genes,
      verbose = verbose
    )
    prior <- network_ref$wk
  }

  # Set initials
  cv_error <- array(0, dim = c(boot_time, 1))
  cv_error_r <- array(0, dim = c(boot_time, 1))
  cv_error_diff <- array(0, dim = c(boot_time, 1))

  # Bootstrap Loop
  for (r in seq_len(boot_time)) {

    # Randomized data
    randomized_data <- as.data.frame(t(apply(data, 1, function(x) {
      sample(x, replace = FALSE)
    })))
    colnames(randomized_data) <- colnames(data)

    # Calculate CV error BIC
    cv_result <- cv_glasso(
      data = data,
      prior = prior,
      k_fold = k_fold,
      crit_cv = crit_cv,
      algorithm = algorithm,
      boot_time = boot_time
    )
    cv_error[r] <- cv_result$avg_error

    cv_result <- cv_glasso(
      data = randomized_data,
      prior = prior,
      k_fold = k_fold,
      crit_cv = crit_cv,
      algorithm = algorithm,
      boot_time = boot_time
    )
    cv_error_r[r] <- cv_result$avg_error
  }

  # CV differential
  cv_error_diff <- (-cv_error) - (-cv_error_r)
  cv <- list(
    cv_error = cv_error, cv_error_r = cv_error_r,
    cv_error_diff = cv_error_diff
  )
  return(cv)
}
