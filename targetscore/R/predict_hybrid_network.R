#' Choose the optimal regulization parameter and scale paramter for prior information adjusted network construction.
#'
#' @param data  input proteomics dataset for network inference(for example: TCGA RPPA data). Gene in coloumns
#' and samples in row. With colnames as gene tags and rownames as sample tags.
#' @param prior Prior information data frame ,with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags. Can be inferred from predict_bio_network() or any network resources.
#' @param cut_off Manually set up cut off value for strength of edge. (Default at 0.1)
#' @param max_dist maximum distance between two antibody. (Default at 1)
#' @param proteomic_responses RPPA data tested for drug pertubation.
#' @param n_prot Antibody number of input data.
#' @param mab_to_genes A list of antibodies, their associated genes, modification sites and effect.
#' @param rho positive tuning parameter vector for elastic net penalty. Default at 10^seq(-2,0, 0.02).
#' @param kappa positive scale parameter vector for prior information matrix contribution. Default at 10^seq(-2,0, 0.02)
#'
#' @return a list is returned with the following entries:
#' {parameters} as the parameter list of regulization parameter decided by the prior information
#' and the algorithmn lowest BIC. Including regularize parameter(L1 norm parameter) as "rho", scale parameter
#' (decided how much prior information contribute) as "kappa", and regulization matrix for the expression
#' data as "rho_m".
#' {bic}{as the Model's BIC error through regularization parameters.}
#' {wk}{inferred network matrix form with edge strength value estimated as the partial correlation.}
#' {wks}{inferred network matrix form with edge strength value estimated as the partial correlation.
#'  Same as wk in predict_hyb_network.}
#' {dist_ind}{A distance file of edgelist with a third column as the network distance between the genes
#'  in the interaction.}
#' {inter}{file as edgelist of inferred network.}
#' {edgelist}{as the edgelist for predicted network.}
#' {nedges}{as the number of edges of the predicted network.}
#'
#' @examples 
#' # READ ANTIBODY FILE ----
#' file <- system.file("targetscoreData", "antibodyMapFile.txt", package = "zeptosensPkg")
#' mab_to_genes <- read.table(file,
#' sep = "\t",
#' header = TRUE,
#' stringsAsFactors = FALSE
#' )
#' 
#' # Read proteomic response for cellline1
#' file <- system.file("test_data", "BT474.csv", package = "zeptosensPkg")
#' proteomic_responses <- read.csv(file, row.names = 1)
#' 
#' # Read Global Signaling file for BRCA
#' file <- system.file("test_data", "TCGA-BRCA-L4.csv", package = "zeptosensPkg")
#' signaling_responses <- read.csv(file, row.names = 1)
#'  
#' # Read Biology knowledge
#' file <- system.file("test_data_files", "predict_bio_network_network_output.rds", 
#'   package = "zeptosensPkg"
#' )
#' prior_org <- readRDS(file)
#'  
#' # Extract network
#' network <- zeptosensPkg::predict_hybrid_network(
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
#' @concept zeptosensPkg
#' @export
predict_hybrid_network <- function(data, prior = NULL, cut_off = 0.1, proteomic_responses, n_prot,
                                   max_dist = 1, mab_to_genes,
                                   rho = 10^seq(-2, 0, 0.02),
                                   kappa = 10^seq(-2, 0, 0.02)) {
  if (is.null(prior)) {
    network_ref <- zeptosensPkg::predict_bio_network(
      n_prot = n_prot,
      proteomic_responses = proteomic_responses,
      max_dist = 1,
      mab_to_genes = mab_to_genes
    )
    wk <- network_ref$wk
    prior <- wk
  }

  # HybNetwork

  index <- colnames(prior[, which(colnames(prior) %in% colnames(data))]) # match the data
  data <- data[, index]
  prior1 <- prior[index, index]

  # prior information extration
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

  # getting the best tuning parameter from BIC minimization
  covmatrix <- (nrow(data) - 1) / nrow(data) * stats::cov(data)
  bic <- matrix(NA, length(rho), length(kappa))
  rho_m <- NULL
  g_result <- NULL
  u <- matrix(1, nrow(prior2), ncol(prior2))
  p_off_d <- NULL
  for (i in seq_len(length(rho))) {
    for (j in seq_len(i)) {
      rho_m <- rho[i] * u - kappa[j] * prior2
      g_result <- glasso::glasso(covmatrix, rho_m, nobs = nrow(covmatrix))
      p_off_d <- sum(g_result$wi != 0 & col(covmatrix) < row(covmatrix))
      bic[i, j] <- -2 * (g_result$loglik) + p_off_d * log(nrow(data))
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
  g_result <- glasso::glasso(covmatrix, rho = rho_m, nobs = nrow(covmatrix))
  sigma_matrix <- g_result$wi
  niter <- g_result$niter
  print(niter) # if niter = 10,000
  if (niter == 10000) {
    stop("ERROR: Algorithmn does not convergence!")
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

  edgelist_total <- zeptosensPkg::create_sif_from_matrix(
    t_net = network_total,
    col_genelist = colnames(network_total),
    row_genelist = rownames(network_total)
  )

  # number of edges
  nedges <- sum(network_total != 0)

  #
  wk <- network_total
  networks <- zeptosensPkg::network2(
    wk = wk, n_prot = n_prot,
    proteomic_responses = proteomic_responses,
  )
  wk <- networks$wk
  wks <- networks$wks
  dist_ind <- networks$dist_ind
  inter <- networks$inter

  # return result
  result <- list(
    parameters = parameters, nedges = nedges, inter = inter,
    wk = wk, wks = wks, dist_ind = dist_ind, edgelist = edgelist_total,
    bic = bic
  )
  return(result)
}
