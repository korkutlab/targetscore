#' Target Score Pilot Run
#'
#' @param wk inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate and -1 for down regulate.
#' @param wks inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate, -1 for down regulate,
#' 2 for phosphorylation and -2 for dephosphorylation.
#' @param dist_ind i distance file of edgelist with a third column as the network distance
#'   between the genes in the interaction
#' @param inter edgelist of inferred network.
#' @param n_dose dose number of input data.
#' @param n_prot antibody number of input data.
#' @param proteomic_responses input drug perturbation data. With columns as antibody, rows as samples.
#' @param max_dist maximum edge strength value.(Default at 1)
#' @param n_perm number of random TS calculations for building the null distribution
#' @param verbose a flag for debugging output
#' @param ts_factor a scaling factor for the pathway component of the target score
#' @param fs_file a file with the functional score data
#'
#' @details
#' data: multiple dose single drug perturbation
#' ts: integral_dose(fs*(xi+sigma_j(2^p*xj*product_k(wk))))
#' to be converted to integral_dose(fs*((xi/(stdev_xi)+sigma_j(2^p*((xj/stdev_xj)*product_k(wk))))
#' missing: For phosp and dephosp based wk, there is no 'exact match' between known and measured phospho-sites
#'
#' @importFrom stats sd pnorm p.adjust
#' @importFrom utils write.table
#'
#' @concept zeptosensPkg
#' @export
get_target_score <- function(wk, wks, dist_ind, inter, n_dose, n_prot, proteomic_responses,
                             max_dist = 1, n_perm, verbose = TRUE, ts_factor = 1,
                             fs_file) {

  # CALCULATE TARGET SCORE ----
  results <- calc_target_score(
    wk = wk,
    wks = wks,
    dist_ind = dist_ind,
    inter = inter,
    n_dose = n_dose,
    n_prot = n_prot,
    proteomic_responses = proteomic_responses,
    max_dist = max_dist,
    verbose = TRUE,
    ts_factor = ts_factor,
    fs_file = fs_file
  )
  ts <- results$ts
  wk <- results$wk
  tsd <- results$tsd
  wks <- results$wks

  # random TS for each node over n permutations comes from randTargetScore.R
  rand_ts <- matrix(0, nrow = n_prot, ncol = n_perm)

  # p value for a given target score computed over the distribution from randTS
  pts <- matrix(0, ncol = 1, nrow = n_prot)

  # CREATE Q-VALUES ----
  for (k in 1:n_perm) {
    #        if(verbose) {
    cat("Permutation Iteration: ", k, "\n")
    #        }

    # print(fs) randomize the readouts over proteomic entities
    rand_proteomic_responses <- proteomic_responses

    # for(j in 1:ncol(rand_proteomic_responses))rand_proteomic_responses[,j] <- sample
    # (proteomic_responses[,j])
    for (i in 1:nrow(rand_proteomic_responses)) rand_proteomic_responses[i, ] <- sample(proteomic_responses[i, ])

    rand_ts[, k] <- calc_target_score(
      wk = wk,
      wks = wks,
      dist_ind = dist_ind,
      inter = inter,
      n_dose = n_dose,
      n_prot = n_prot,
      proteomic_responses = rand_proteomic_responses,
      max_dist = max_dist,
      cell_line = cell_line,
      verbose = verbose,
      ts_factor = ts_factor,
      fs_file = fs_file
    )$ts

    # rand_ts[,k] <- as.matrix(rants) print('resi') print(resi$ts) rand_ts[,k]
  }

  for (i in 1:n_prot) {
    mean <- mean(rand_ts[i, 1:n_perm])
    stdev <- sd(rand_ts[i, 1:n_perm])
    zval <- (ts[i] - mean) / (stdev)
    pts[i] <- 2 * pnorm(-abs(zval)) # pnorm(ts[i], mean = mean(rand_ts[i, 1:n_perm]), sd = sd(rand_ts[i, 1:n_perm]))

    if (verbose) {
      print(pts[i])
    }
  }

  q <- as.matrix(p.adjust(pts, method = "fdr", n = n_prot))

  # Set row and column names for results
  rownames(q) <- colnames(proteomic_responses)
  colnames(q) <- "FDR_adjusted_p"

  rownames(rand_ts) <- colnames(proteomic_responses)
  rownames(pts) <- colnames(proteomic_responses)

  # RETURN RESULTS ----
  results <- list(ts = ts, wk = wk, tsd = tsd, q = q, wks = wks, pts = pts, rand_ts = rand_ts)

  return(results)
}
