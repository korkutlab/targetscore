#' Choose the optimal regulization parameter and scale paramter for adjusted-glasso algorithm network construction
#'
#' @param data  input proteomics dataset for network inference. Gene in coloumns and samples in row.
#' With colnames as gene tags and rownames as sample tags.
#' @param prior Prior information matrix of gene interaction, with colnames and rownames as gene tags.
#' With colnames and rownames as gene tags.Can be inferred from Public data source (for example:SignedPC).
#' @param rho positive tuning parameter for elastic net penalty. Default at seq10^seq(-2, 0, 0.02).
#' @param kappa positive scaler parameter for biology-knowledge base contribution. Default at seq10^seq(-2, 0, 0.02).
#'
#' @return Parameter list of regulization parameter decided by the prior information and the algorithmn lowest BIC.
#' \item{rho}{optimized penalty parameter.}
#' \item{kappa}{optimized scaler parameter.}
#' \item{rho_m}{optimized regulization matrix for the expression data.}
#' \item{bic}{Calculated model BIC error for differnet regularization parameters.}
#' 
#' @examples 
#' # read proteomic response file
#' file <- system.file("test_data", "TCGA-BRCA-L4.csv", package = "zeptosensPkg")
#' signaling_responses <- read.csv(file, row.names = 1)
#' 
#' # Read in Biology knowlegde base protein interaction
#' file <- system.file("test_data_files", "predict_bio_network_network_output.rds",
#'   package = "zeptosensPkg")
#' 
#' prior_org <- readRDS(file)
#' parameters <- zeptosensPkg::optimize_parameter_hybrid(data = signaling_responses, 
#'   prior = prior_org$wk)
#' 
#' @importFrom glasso glasso
#' @importFrom stats cov
#'
#' @concept zeptosensPkg
#' @export
optimize_parameter_hybrid <- function(data, prior = NULL,
                                      rho = 10^seq(-2, 0, 0.02),
                                      kappa = 10^seq(-2, 0, 0.02)) {
  # Extract from SignedPC for prior

  # READ ANTIBODY FILE ----
  mab_to_genes <- read.table(system.file("targetscoreData", "antibodyMapFile_08092019.txt", package = "zeptosensPkg"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )

  if (is.null(prior)) {
    network_ref <- zeptosensPkg::predict_bio_network(
      n_prot = ncol(data),
      proteomic_responses = data,
      max_dist = 1,
      mab_to_genes = mab_to_genes
    )
    wk <- network_ref$wk
    prior <- wk
  }

  # Match the data with prior
  index <- colnames(prior[, which(colnames(prior) %in% colnames(data))]) # match the data
  data <- data[, index]
  prior1 <- prior[index, index]

  # Symetric prior information extration
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

  # Calculate covariance matrix
  covmatrix <- (nrow(data) - 1) / nrow(data) * stats::cov(data)

  # Set initial value
  bic <- matrix(NA, length(rho), length(kappa))
  rho_m <- NULL
  g_result <- NULL
  u <- matrix(1, nrow(prior2), ncol(prior2))
  p_off_d <- NULL

  # Getting the best tuning parameter from BIC minimization

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
  parameters <- list(rho_m = rho_m, rho = rho, kappa = kappa, bic = bic)


  return(parameters)
}
