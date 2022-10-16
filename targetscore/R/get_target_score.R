#' Get Target Score Value Significance
#'
#' @param wk inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate and -1 for down regulate.
#' @param wks inference network constructed in matrix form with edge strength value estimated.
#' Can be extracted directly from predict_bio_network,predict_dat_network or predictt_hyb_network
#' function. Where predict_bio_network edge value default at 1 for upregulate, -1 for down regulate,
#' 2 for phosphorylation and -2 for dephosphorylation.
#' @param dist_ind i distance file of edgelist with a third column as the network distance
#' between the genes in the interaction
#' @param edgelist edgelist of inferred network.
#' @param n_dose dose number of input data.
#' @param n_prot antibody number of input data.
#' @param proteomic_responses input drug perturbation data. With columns as antibody, rows as samples.
#' @param n_perm number of random TS calculations for building the null distribution
#' @param verbose a flag for debugging output
#' @param ts_pathway_scale a scaling factor for the pathway component of the target score
#' @param fs_dat a dataset with the functional score data.First coloumn as the protein name and second
#' column as the functional score. Can be inferred from get_fs_value or be User defined and curated.
#' 
#' @return a list is returned with the following entries:
#' {ts}{TargetScore values summed over individual drug doses}
#' {tsd}{TargetScore values for individual drug doses}
#' {wk}{Inferred network matrix form with edge strength value estimated as the partial correlation}
#' {wks}{Inferred network matrix form with edge strength value estimated as the partial correlation}
#' {pts}{Unadjusted p-values calculated for using permutation testing where proteomic responses for each row have been randomized}
#' {q}{Q-values (derived from FDR adjusted p-values) calculated for using permutation testing where proteomic responses for each row have been randomized}
#' {rand_ts}{TargetScores for each randomized permutation}
#' 
#' @examples 
#' # read proteomic response file
#' signaling_responses <- read.csv(system.file("test_data", "TCGA-BRCA-L4.csv",
#' package = "targetscore"
#' ), row.names = 1)
#' 
#' # Read in biology knowledge base protein interaction
#' network <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
#' package = "targetscore"
#' ))
#' 
#' # read proteomic responce file
#' proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"),
#' row.names = 1
#' )
#' 
#' # read functional score value
#' fs <- readRDS(system.file("test_data_files", "get_fs_vals_output.rds",
#' package = "targetscore"
#' ))
#' 
#' # Calculate Target Score
#' ts <- array(0, dim = c(dim(proteomic_responses)[1], dim(proteomic_responses)[2]))
#' ts_p <- array(0, dim = c(dim(proteomic_responses)[1], dim(proteomic_responses)[2]))
#' ts_q <- array(0, dim = c(dim(proteomic_responses)[1], dim(proteomic_responses)[2]))
#' for (i in seq_len(2)) {
#' results <- targetscore::get_target_score(
#'   wk = network$wk,
#'   wks = network$wks,
#'   dist_ind = network$dist_ind,
#'   edgelist = network$edgelist,
#'   n_dose = 1,
#'   n_prot = dim(proteomic_responses)[2],
#'   proteomic_responses = proteomic_responses[i,],
#'   n_perm = 1,
#'   verbose = FALSE,
#'   fs_dat = fs
#' )
#' ts[i,] <- t(results$ts)
#' ts_p[i,] <- t(results$pts)
#' ts_q[i,] <- t(results$q)
#' } 
#' colnames(ts) <- colnames(proteomic_responses)
#' ts <- data.frame(rownames(proteomic_responses), ts)
#' colnames(ts_p) <- colnames(proteomic_responses)
#' ts_p <- data.frame(rownames(proteomic_responses), ts_p)
#' colnames(ts_q) <- colnames(proteomic_responses)
#' ts_q <- data.frame(rownames(proteomic_responses), ts_q)
#' ts_result <- list(ts = ts, ts_p = ts_p, ts_q = ts_q)
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
#' @concept targetscore
#' @export
get_target_score <- function(wk, wks, dist_ind, edgelist, n_dose, n_prot, proteomic_responses,
                             n_perm, verbose = TRUE, ts_pathway_scale = 1, fs_dat) {

  # CALCULATE TARGET SCORE ----
  results <- calc_target_score(
    wk = wk,
    wks = wks,
    dist_ind = dist_ind,
    edgelist = edgelist,
    n_dose = n_dose,
    n_prot = n_prot,
    proteomic_responses = proteomic_responses,
    verbose = verbose,
    ts_pathway_scale = ts_pathway_scale,
    fs_dat = fs_dat
  )
  
  ts <- results$ts
  wk <- results$wk
  tsd <- results$tsd
  wks <- results$wks

  # Random TS for each node over n permutations
  rand_ts <- matrix(NA, nrow = n_prot, ncol = n_perm)

  # p value for a given target score computed over the distribution from randTS
  pts <- matrix(NA, ncol = 1, nrow = n_prot)
  
  ## This will be the final targetscore summed over multiple doses; always 1 row
  #ts <- matrix(NA, ncol = n_prot, nrow = 1)

  # CREATE Q-VALUES ----
  for (k in seq_len(n_perm)) {
    #        if(verbose) {
    message("MSG: Permutation Iteration: ", k, "\n")
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
      edgelist = edgelist,
      n_dose = n_dose,
      n_prot = n_prot,
      proteomic_responses = rand_proteomic_responses,
      verbose = verbose,
      ts_pathway_scale = ts_pathway_scale,
      fs_dat = fs_dat
    )$ts

    # rand_ts[,k] <- as.matrix(rants) print('resi') print(resi$ts) rand_ts[,k]
  }

  # NaN are created because the ts, mean, and sd are 0
  for (i in 1:n_prot) {
    mean <- mean(rand_ts[i, 1:n_perm])
    stdev <- sd(rand_ts[i, 1:n_perm])
    zval <- (ts[i] - mean) / (stdev)
    pts[i] <- 2 * pnorm(-abs(zval)) # pnorm(ts[i], mean = mean(rand_ts[i, 1:n_perm]), sd = sd(rand_ts[i, 1:n_perm]))

    if (verbose) {
      print(pts[i])
    }
  }

  # q and pts can be named vectors; ts must be a data.frame for multiple entries
  q <- as.matrix(p.adjust(pts, method = "fdr", n = n_prot))
  q <- t(q)

  # Set row and column names for results
  colnames(q) <- colnames(proteomic_responses)

  # NOTE: At this point in the code all the outputs have been generated
  # rand_ts: a matrix; (TS number by permutation count) 
  # pts: a vector; 1 TS for each entity
  # q: a vector 1 TS for each entity
  # ts: a vector 1 TS for each entity
  
  # Process outputs
  rand_ts <- t(rand_ts)
  colnames(rand_ts) <- colnames(proteomic_responses)
  
  pts <- t(pts)
  colnames(pts) <- colnames(proteomic_responses)
  
  ts <- data.frame(as.list(ts))
  colnames(ts) <- colnames(proteomic_responses)
  row.names(ts) <- NULL
  
  # tsd is in the correct format
  
  # # Make sure everything is a named vector before returning results
  # ts <- data.frame(as.list(ts))
  # colnames(ts) <- colnames(proteomic_responses)
  # ts <- unlist(ts)
  # 
  # tsd <- data.frame(as.list(tsd))
  # colnames(tsd) <- colnames(proteomic_responses)
  # tsd <- unlist(tsd)
  # 
  # pts <- data.frame(as.list(pts))
  # colnames(pts) <- colnames(proteomic_responses)
  # pts <- unlist(pts)
  # 
  #q_vec <- data.frame(as.list(q))
  #colnames(q_vec) <- colnames(proteomic_responses)
  #q_vec <- unlist(q)
  
  # NOTE: This will not work because ts can be multiple columns
  ts_sig_df <- t(rbind(ts, pts, q))
  colnames(ts_sig_df) <- c("ts", "p_value", "q_value")
  ts_sig_df <- as.data.frame(ts_sig_df) # Returns a matrix with without as.data.frame
  
  # RETURN RESULTS ----
  results <- list(wk=wk, wks=wks, 
                  ts=ts, tsd=tsd, rand_ts=rand_ts,
                  pts=pts, q=q, ts_sig_df=ts_sig_df)
  
  return(results)
}
