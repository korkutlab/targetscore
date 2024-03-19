#' Choose the optimal regularization parameter and scale parameter for prior information adjusted network construction.
#'
#' @param data input proteomics dataset for network inference (for example: TCGA RPPA data). Gene in columns
#' and samples in row. With colnames as gene tags and rownames as sample tags.
#' @param prior prior information data frame, with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags. Can be inferred from predict_bio_network() or any network resources.
#' @param cut_off manually set up cut off value for strength of edge. (Default at 0.1)
#' @param max_dist maximum distance between two antibody. (Default at 1)
#' @param proteomic_responses Input proteomics dataset for network inference (for 
#'   example: TCGA RPPA data or post-perturbation response data).
#'   With columns as antibody, rows as samples (Defaults to "data")
#' @param n_prot antibody number of input data.
#' @param mab_to_genes a list of antibodies, their associated genes, modification sites and effect.
#' @param rho positive tuning parameter vector for elastic net penalty. Default at 10^seq(-2, 0, 0.02).
#' @param kappa positive scale parameter vector for prior information matrix contribution. Default at 10^seq(-2, 0, 0.02)
#' @param verbose logical, whether to show additional debugging information
#' 
#' @note proteomic_responses is used only to retrieve the desired list of 
#' entries for the resulting network
#'
#' @return a list is returned with the following entries:
#' {parameters} as the parameter list of regularization parameter decided by the prior information
#'   and the algorithm lowest BIC. Including regularize parameter(L1 norm parameter) as "rho", scale parameter
#'   (decided how much prior information contribute) as "kappa", and regulization matrix for the expression
#'   data as "rho_m".
#' {bic}{as the model's BIC error through regularization parameters.}
#' {wk}{inferred network matrix form with edge strength value estimated as the partial correlation.}
#' {wks}{inferred network matrix form with edge strength value estimated as the partial correlation.
#'   Same as wk in predict_hyb_network.}
#' {dist_ind}{A distance file of edgelist with a third column as the network distance between the genes
#'   in the interaction.}
#' {inter}{file as edgelist of inferred network.}
#' {edgelist}{as the edgelist for predicted network.}
#' {nedges}{as the number of edges of the predicted network.}
#'
#' @examples 
#' # READ ANTIBODY FILE ----
#' file <- system.file("target_score_data", "antibody_map.csv", package = "targetscore")
#' mab_to_genes <- read.csv(file,
#' header = TRUE,
#' stringsAsFactors = FALSE
#' )
#' 
#' # Read proteomic response for cellline1
#' file <- system.file("test_data", "BT474.csv", package = "targetscore")
#' proteomic_responses <- read.csv(file, row.names = 1)
#' 
#' # Read Global Signaling file for BRCA
#' file <- system.file("test_data", "TCGA-BRCA-L4.csv", package = "targetscore")
#' signaling_responses <- read.csv(file, row.names = 1)
#'  
#' # Read Biology knowledge
#' file <- system.file("test_data_files", "predict_bio_network_network_output.rds", 
#'   package = "targetscore"
#' )
#' prior_org <- readRDS(file)
#'  
#' # Extract network
#' network <- targetscore::predict_hybrid_network(
#' data = signaling_responses,
#' prior = prior_org$wk,
#' n_prot = dim(proteomic_responses)[2],
#' proteomic_responses = proteomic_responses,
#' mab_to_genes = mab_to_genes,
#' max_dist = 1
#' )
#'
#' @importFrom glasso glasso
#' @importFrom stats cov median na.omit
#'
#' @concept targetscore
#' @export
predict_hybrid_network <- function(data, prior = NULL, cut_off = 0.1, 
                                   proteomic_responses=data, n_prot,
                                   max_dist = 1, mab_to_genes,
                                   rho = 10^seq(-2, 0, 0.02),
                                   kappa = 10^seq(-2, 0, 0.02),
                                   verbose=FALSE) {
  if (is.null(prior)) {
    network_ref <- targetscore::predict_bio_network(
      n_prot = n_prot,
      proteomic_responses = proteomic_responses,
      max_dist = 1,
      mab_to_genes = mab_to_genes
    )
    wk <- network_ref$wk
    prior <- wk
  }

  # Hybrid Network
  index <- colnames(prior[, which(colnames(prior) %in% colnames(data))]) # match the data
  data <- data[, index]
  prior1 <- prior[index, index]

  # Prior information extraction
  prior1 <- ifelse(prior1 != 0, 1, 0) # information matrix of prior
  prior2 <- prior1 # symmetrical prior information
  for (i in seq_len(nrow(prior1))) {
    for (j in seq_len(ncol(prior1))) {
      if (prior1[i, j] != 0) {
        prior2[i, j] <- prior1[i, j]
        prior2[j, i] <- prior1[i, j]
      }
    }
  }
  prior2 <- ifelse(prior2 != 0, 1, 0)

  # FIXME: FIG 3B BIC MATRIX 
  # getting the best tuning parameter from BIC minimization
  covmatrix <- (nrow(data) - 1) / nrow(data) * stats::cov(data)
  bic <- matrix(NA, length(rho), length(kappa))
  rho_m <- NULL
  g_result <- NULL
  u <- matrix(1, nrow(prior2), ncol(prior2))
  p_off_d <- NULL # FIXME: DIAGONAL ELEMENTS
  for (i in seq_len(length(rho))) {
    for (j in seq_len(i)) {
      rho_m <- rho[i] * u - kappa[j] * prior2
      g_result <- glasso::glasso(covmatrix, rho_m, nobs = nrow(covmatrix))
      p_off_d <- sum(g_result$wi != 0 & col(covmatrix) < row(covmatrix))
      bic[i, j] <- -2 * (g_result$loglik) + p_off_d * log(nrow(data)) # FIXME: BIC ERROR?; see: https://www.stata.com/meeting/us21/slides/US21_Dallakyan.pdf
      bic <- as.data.frame(bic)
      rownames(bic) <- rho
      colnames(bic) <- kappa
    }
  }
  
  pos <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)
  rho <- rho[pos[1]]
  kappa <- kappa[pos[2]]
  rho_m <- rho * u - kappa * prior2
  parameters <- list(rho_m = rho_m, rho = rho, kappa = kappa)

  # Estimated inverse covariance (precision)
  # Default glasso max iterations (maxit) is 10000
  g_result <- glasso::glasso(covmatrix, rho = rho_m, nobs = nrow(covmatrix))
  sigma_matrix <- g_result$wi
  niter <- g_result$niter
  
  if(verbose) {
    print(niter) # if niter = 10,000    
  }
  
  if (niter == 10000) {
    stop("ERROR: Algorithm does not converge!")
  }

  pcor_matrix <- matrix(0, nrow = ncol(data), ncol = ncol(data))
  
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      pcor_matrix[i, j] <- -sigma_matrix[i, j] / sqrt(sigma_matrix[i, i] * sigma_matrix[j, j])
    }
  }

  t_edges <- pcor_matrix

  # cut off = cut off value
  t_net <- as.data.frame(ifelse(abs(t_edges) >= cut_off & row(t_edges) != col(t_edges), t_edges, 0))
  colnames(t_net) <- colnames(data)
  rownames(t_net) <- colnames(data)

  # directed network
  # get direction for the network as the bigger covariance estimated indicated the upper stream gene;
  t_net_d <- matrix(0, nrow = nrow(t_net), ncol = ncol(t_net))
  for (i in 1:ncol(t_net)) {
    for (j in 1:nrow(t_net)) {
      if (prior1[i, j] != 0 & prior1[j, i] == 0) {
        t_net_d[i, j] <- t_net_d[i, j]
      }
      if (prior1[i, j] == 0 & prior1[j, i] == 0) {
        t_net_d[i, j] <- t_net[i, j]
        t_net_d[j, i] <- t_net[j, i]
      }
      if (prior1[i, j] != 0 & prior1[j, i] != 0) {
        t_net_d[i, j] <- t_net[i, j]
        t_net_d[j, i] <- t_net[j, i]
      }
    }
  }
  colnames(t_net_d) <- colnames(data)
  rownames(t_net_d) <- colnames(data)

  # Inferred missing proteins from bio-network
  if (ncol(prior) > ncol(prior1)) {
    network <- t_net_d
    index2 <- which(colnames(prior) %in% index)
    prior1_extra <- prior[-index2, -index2]
    index_extra <- colnames(prior1_extra)

    # Adding the missing node from prior with the edgevalue
    ## (prior information have the value of the median(abd(netowrk)))
    prior3 <- as.data.frame(prior * median(abs(c(na.omit(c(ifelse(network != 0, network, NA)))))))
    colnames(prior3) <- colnames(prior)
    rownames(prior3) <- rownames(prior)

    part1 <- prior3[which(rownames(prior3) %in% index_extra), which(colnames(prior) %in% index)]
    part2 <- prior3[which(rownames(prior3) %in% index), colnames(prior) %in% index_extra]
    part3 <- prior3[which(rownames(prior3) %in% index_extra), which(colnames(prior) %in% index_extra)]

    network_total <- cbind(rbind(network, part1), rbind(part2, part3))
    index3 <- match(colnames(prior), colnames(network_total))
    network_total <- network_total[index3, index3]
    colnames(network_total) <- colnames(prior)
    rownames(network_total) <- rownames(prior)
  }

  if (length(prior) == length(prior1)) {
    network_total <- t_net_d
  }

  edgelist_total <- targetscore::create_sif_from_matrix(
    t_net = network_total,
    col_genelist = colnames(network_total),
    row_genelist = rownames(network_total)
  )

  # Number of edges
  nedges <- sum(network_total != 0)

  wk <- network_total
  networks <- targetscore::predict_dat_network_get_properties(
    wk = wk, n_prot = n_prot,
    proteomic_responses = proteomic_responses,
  )
  wk <- networks$wk
  wks <- networks$wks
  dist_ind <- networks$dist_ind
  inter <- networks$inter

  # Return result
  result <- list(
    parameters = parameters, nedges = nedges, inter = inter,
    wk = wk, wks = wks, dist_ind = dist_ind, edgelist = edgelist_total,
    bic = bic
  )
  
  return(result)
}
