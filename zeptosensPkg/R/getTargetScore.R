#' Target Score Pilot Run
#'
#' @param wk TODO
#' @param wks TODO
#' @param dist_ind TODO
#' @param inter TODO
#' @param n_dose TODO
#' @param n_prot TODO
#' @param proteomic_responses TODO
#' @param max_dist TODO (default: 1)
#' @param cell_line TODO
#' @param target_score_output_file a filename to write total target score results (default: NULL)
#' @param matrix_wk_output_file TODO
#' @param target_score_q_value_file a filename to write statistical significance levels (default: NULL)
#' @param target_score_p_value_file a filename to write statistical significance levels (default: NULL)
#' @param target_score_dose_file a filename to write dose dependent target score results (default: NULL)
#' @param n_perm number of random TS calculations for building the null distribution
#' @param verbose a flag for debugging output
#' @param ts_factor a scaling factor for the pathway component of the target score
#' @param fs_file a file with the functional score data
#' @param signed_matrix_wk_output_file TODO
#' @param random_target_score_file TODO
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
getTargetScore <- function(wk, wks, dist_ind, inter, n_dose, n_prot, proteomic_responses,
                           max_dist = 1, n_perm, cell_line, target_score_output_file = NULL,
                           matrix_wk_output_file = NULL, target_score_q_value_file = NULL,
                           target_score_dose_file = NULL, random_target_score_file = NULL,
                           target_score_p_value_file = NULL, verbose = TRUE, ts_factor = 1,
                           fs_file, signed_matrix_wk_output_file = NULL) {

  # CALCULATE TARGET SCORE ----
  results <- calcTargetScore(
    wk = wk,
    wks = wks,
    dist_ind = dist_ind,
    inter = inter,
    n_dose = n_dose,
    n_prot = n_prot,
    proteomic_responses = proteomic_responses,
    max_dist = max_dist,
    cell_line = cell_line,
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

    rand_ts[, k] <- calcTargetScore(
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
  rownames(q) <- colnames(proteomic_responses)
  colnames(q) <- "FDR_adjusted_p"

  # WRITE OUTPUTS ----
  if (!is.null(matrix_wk_output_file)) {
    write.table(wk, file = matrix_wk_output_file, quote = FALSE, sep = "\t")
  }
  if (!is.null(signed_matrix_wk_output_file)) {
    write.table(wks, file = signed_matrix_wk_output_file, quote = FALSE, sep = "\t")
  }
  if (!is.null(target_score_output_file)) {
    write.table(ts, file = target_score_output_file, quote = FALSE, col.names = FALSE, sep = "\t")
  }

  if (!is.null(target_score_dose_file)) {
    write.table(tsd, file = target_score_dose_file, quote = FALSE, sep = "\t")
  }

  if (!is.null(random_target_score_file)) {
    write.table(data.frame(rand_ts), file = random_target_score_file, quote = FALSE, sep = "\t")
  }

  if (!is.null(target_score_q_value_file)) {
    write.table(q, file = target_score_q_value_file, quote = FALSE, sep = "\t")
  }
  rownames(rand_ts) <- colnames(proteomic_responses)
  rownames(pts) <- colnames(proteomic_responses)

  write.table(rand_ts, file = "randts.txt", quote = FALSE, sep = "\t")
  write.table(pts, file = target_score_p_value_file, quote = FALSE, sep = "\t")

  # RETURN RESULTS ----
  results <- list(ts = ts, wk = wk, tsd = tsd, q = q, wks = wks, pts = pts)

  return(results)
}
