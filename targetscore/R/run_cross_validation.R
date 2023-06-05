#' Run cross validation through data.
#'
#' @param data input proteomics dataset for network inference. Gene in columns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param prior prior information matrix of gene interaction, with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags. Can be inferred from Public data source (for example:SignedPC).
#' @param n_prot antibody number of input data.
#' @param k_fold specify the number of folds for cross validation. (Default at 5)
#' @param proteomic_responses input drug perturbation data. With columns as antibody, rows as samples.
#' @param mab_to_genes A list of antibodies, their associated genes, modification sites and effect.
#' @param dist_file a distance file an edgelist with a third column which is the network distance
#' between the genes in the interaction
#' @param max_dist maximum distance between two antibody. (Default at 1)
#' @param verbose whether to show debugging information
#' @param boot_time bootstrap time manually set. Default at 100.
#' @param crit_cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param algorithm flexible toolbox implementing network estimating algorithms for robustness test.
#' (\code{data_driven},or \code{hybrid_driven}).
#'
#' @return returns list of cross validation results which includes:
#' \item{cv_error}{cross validation errors of data.}
#' \item{cv_error_r}{cross validation errors of random data.}
#' \item{cv_error_diff}{cross validation errors difference between random data and data.}
#'
#' @importFrom glasso glasso
#'
#' @concept targetscore
#' @export
#' 
#' @details FIXME: FIG 3D BIORXIV
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
                                 crit_cv = c("loglik", "AIC", "BIC"),
                                 algorithm = c("data_driven", "hybrid_driven")) {
  
  # BASIC CHECK ----
  if(!all(colnames(data) %in% mab_to_genes$AntibodyLabel)) {
    stop("ERROR: Not all 'data' columns are in antibody map")
  }
  
  if(!all(colnames(proteomic_responses) %in% mab_to_genes$AntibodyLabel)) {
    stop("ERROR: Not all 'proteomic_responses' columns are in antibody map")  
  }
  
  ## Test for missing data
  missing_val_cnt <- apply(data, 2, function(x) { length(which(is.na(x))) }) 
  if(any(missing_val_cnt != 0)) {
    stop("ERROR: Some 'data' columns have missing data")  
  }
  
  # MATCH DATA ----
  # Labels need to match in "data" and "proteomic_responses"
  #  "proteomic_responses" should use labels from "data" 
  
  old_cols <- proteomic_responses
  new_cols <- data
  
  new_colnames <- rep(NA, ncol(old_cols))
  
  for(i in 1:ncol(old_cols)) {
    #i <- 1
    #cat("I: ", i, "\n")
    cur_colname <- colnames(old_cols)[i]
    cur_col_idx <- which(mab_to_genes$AntibodyLabel == cur_colname)
    tmp_col <- mab_to_genes[cur_col_idx,]
    
    tmp_mab_to_genes <- mab_to_genes[which(mab_to_genes$AntibodyLabel %in% colnames(new_cols)),]
    
    b <- which(tmp_mab_to_genes$Gene_Symbol %in% tmp_col$Gene_Symbol)
    if(is.na(unique(tmp_col$Sites))) {
      c <- which(is.na(tmp_mab_to_genes$Sites))
    } else {
      c <- which(tmp_mab_to_genes$Sites == unique(tmp_col$Sites))
    }
    d <- which(tmp_mab_to_genes$Effect == unique(tmp_col$Effect))
    idx <- Reduce(intersect, list(b, c, d))
    
    antibody_label <- unique(tmp_mab_to_genes$AntibodyLabel[idx])
    
    if(length(idx) == 1) {
      new_colnames[i] <- antibody_label
    } 
  }
  
  ncol(old_cols)
  idx <- which(!is.na((new_colnames)))
  new_proteomics_data <- old_cols[,idx]
  colnames(new_proteomics_data) <- new_colnames[idx]
  
  # head(new_proteomics_data)
  # head(data)
  
  proteomics_data <- new_proteomics_data

  # Old code
  # index <- colnames(data[which(colnames(data) %in% colnames(proteomic_responses))])
  # data[, index]
  # proteomic_responses <- proteomic_responses[, index]

  # Import prior biology information/NULL from SignedPC
  if (is.null(prior)) {
    network_ref <- targetscore::predict_bio_network(
      n_prot = n_prot,
      proteomic_responses = proteomic_responses,
      max_dist = max_dist,
      mab_to_genes = mab_to_genes,
      verbose = verbose
    )
    
    prior <- network_ref$wk
  }
  
  # DEBUG
  # NOTE: To debug errors that say: arguments imply differing number of rows: 0, 1
  # setdiff(colnames(proteomic_responses), mab_to_genes$AntibodyLabel)

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
    )
    cv_error[r] <- cv_result$avg_error

    cv_result <- cv_glasso(
      data = randomized_data,
      prior = prior,
      k_fold = k_fold,
      crit_cv = crit_cv,
      algorithm = algorithm,
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
