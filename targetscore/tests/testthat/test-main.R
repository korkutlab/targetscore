test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("match_genes_to_edgelist", {
  antibody_map_file <- system.file("target_score_data", "antibody_map.csv", package = "targetscore")
  mab_to_genes <- read.csv(antibody_map_file, header = TRUE, stringsAsFactors = FALSE)
  
  proteomic_responses_file <- system.file("test_data", "BT474.csv", package = "targetscore")
  proteomic_responses <- read.csv(proteomic_responses_file, row.names = 1)

  dist_file <- system.file("target_score_data", "distances.txt", package = "targetscore")
  tmp_dist <- read.table(dist_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  dist <- tmp_dist[which(tmp_dist[, 3] <= 1), ]

  idx_ab_map <- which(mab_to_genes[, 1] %in% colnames(proteomic_responses))

  mab_genes <- mab_to_genes[idx_ab_map, 4]
  names(mab_genes) <- mab_to_genes[idx_ab_map, 1]

  dist_list <- match_genes_to_edgelist(
    genes1 = mab_genes,
    genes2 = NULL,
    annot_edgelist = dist,
    antibody_vec = colnames(proteomic_responses),
    use_annot = TRUE,
    verbose = TRUE
  )

  snapshot_file <- system.file("test_data_files", "match_genes_to_edgelist_test_output.rds", package = "targetscore")
  #saveRDS(dist_list, "match_genes_to_edgelist_test_output.rds")
  snapshot <- readRDS(snapshot_file)

  expect_identical(dist_list, snapshot)
})

test_that("predict_bio_network", {
  skip_on_cran()

  # READ ANTIBODY FILE ----
  mab_to_genes <- read.csv(system.file("target_score_data", "antibody_map.csv", package = "targetscore"),
    header = TRUE,
    stringsAsFactors = FALSE
  )

  # Read proteomic response for cellline1
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"), row.names = 1)

  # Extract network
  network <- targetscore::predict_bio_network(
    n_prot = dim(proteomic_responses)[2],
    proteomic_responses = proteomic_responses,
    max_dist = 1,
    mab_to_genes = mab_to_genes
  )

  network_org <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
    package = "targetscore"
  ))

  expect_identical(network$wk, network_org$wk)
  expect_identical(network$wks, network_org$wks)
  expect_identical(network$dist_ind, network_org$dist_ind)
  expect_identical(network$inter, network_org$inter)
})

test_that("predict_dat_network", {
  skip_on_cran()

  # Read proteomic response for cellline1
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"), row.names = 1)

  # Read Global Signaling file for BRCA
  signaling_responses <- read.csv(system.file("test_data", "TCGA-BRCA-L4.csv", package = "targetscore"), row.names = 1)

  # Extract network
  network <- targetscore::predict_dat_network(
    data <- signaling_responses,
    n_prot = dim(proteomic_responses)[2],
    proteomic_responses = proteomic_responses
  )

  network_org <- readRDS(system.file("test_data_files", "predict_dat_network_output.rds",
    package = "targetscore"
  ))

  expect_identical(network$wk, network_org$wk)
  expect_identical(network$wks, network_org$wks)
  expect_identical(network$dist_ind, network_org$dist_ind)
  expect_identical(network$inter, network_org$inter)
})

test_that("predict_hybrid_network", {
  skip_on_cran()

  # READ ANTIBODY FILE ----
  mab_to_genes <- read.csv(system.file("target_score_data", "antibody_map.csv", package = "targetscore"),
    header = TRUE,
    stringsAsFactors = FALSE
  )

  # Read proteomic response for cellline1
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"), row.names = 1)

  # Read Global Signaling file for BRCA
  signaling_responses <- read.csv(system.file("test_data", "TCGA-BRCA-L4.csv", package = "targetscore"), row.names = 1)

  # Read Biology knowledge
  prior_org <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
    package = "targetscore"
  ))

  # Extract network
  network <- targetscore::predict_hybrid_network(
    data = signaling_responses,
    prior = prior_org$wk,
    n_prot = dim(proteomic_responses)[2],
    proteomic_responses = proteomic_responses,
    mab_to_genes = mab_to_genes,
    max_dist = 1
  )

  network_org <- readRDS(system.file("test_data_files", "predict_hybrid_network_network_output.rds",
    package = "targetscore"
  ))

  expect_identical(network$wk, network_org$wk)
  expect_identical(network$wks, network_org$wks)
  expect_identical(network$dist_ind, network_org$dist_ind)
  expect_identical(network$inter, network_org$inter)
})

test_that("predict_dat_network_get_properties", {

  # Read proteomic response for cellline1
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"), row.names = 1)

  # Read network output
  network_org <- readRDS(system.file("test_data_files", "network_properties_output.rds",
    package = "targetscore"
  ))

  # Read in network output
  wk_org <- readRDS(system.file("test_data_files", "predict_hybrid_network_network_output.rds",
    package = "targetscore"
  ))

  network <- targetscore::predict_dat_network_calc_properties(
    wk <- wk_org$wk,
    n_prot = dim(proteomic_responses)[2],
    proteomic_responses = proteomic_responses
  )

  expect_identical(network$wk, network_org$wk)
  expect_identical(network$wks, network_org$wks)
  expect_identical(network$dist_ind, network_org$dist_ind)
  expect_identical(network$inter, network_org$inter)
})


test_that("create_sif_from_matrix", {
  network_org <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
    package = "targetscore"
  ))

  wk_org <- readRDS(system.file("test_data_files", "create_sif_from_matrix_output.rds",
    package = "targetscore"
  ))

  edgelist_wk <- create_sif_from_matrix(
    t_net = network_org$wk,
    col_genelist = colnames(network_org$wk),
    row_genelist = rownames(network_org$wk)
  )
  
  # Save snapshot
  #saveRDS(edgelist_wk, file.path("inst", "test_data_files", "create_sif_from_matrix_output.rds"))

  expect_identical(wk_org, edgelist_wk)
})

test_that("get_fs_vals", {
  # Read fs_manually set file
  fs_override_org <- readRDS(system.file("test_data_files", "fs_value_file.rds",
    package = "targetscore"
  ))

  # Read proteomic response file
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"), row.names = 1)
  
  # Read antibody file
  mab_to_genes <- read.csv(system.file("target_score_data", "antibody_map.csv", package = "targetscore"),
    header = TRUE,
    stringsAsFactors = FALSE
  )

  wk_org <- readRDS(system.file("test_data_files", "get_fs_vals_output.rds",
    package = "targetscore"
  ))

  fs <- get_fs_vals(
    n_prot = ncol(proteomic_responses), proteomic_responses = proteomic_responses,
    mab_to_genes = mab_to_genes, fs_override = fs_override_org
  )

  expect_identical(fs, wk_org)
})

test_that("samp_sdev", {

  # read proteomic responce file
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"), row.names = 1)

  # test-data sdev
  sd_org <- readRDS(system.file("test_data_files", "samp_sdev_output.rds",
    package = "targetscore"
  ))

  samp_d <- samp_sdev(
    n_x = proteomic_responses, n_sample = dim(proteomic_responses)[1],
    n_prot = dim(proteomic_responses)[2], n_dose = 1
  )

  expect_identical(samp_d, sd_org)
})

test_that("optimize_parameter_dat", {

  # read proteomic response file
  signaling_responses <- read.csv(system.file("test_data", "TCGA-BRCA-L4.csv", package = "targetscore"), row.names = 1)

  # test-data data parameter output
  parameter_org <- readRDS(system.file("test_data_files", "optimize_parameter_dat_output.rds",
    package = "targetscore"
  ))
  parameters <- targetscore::optimize_parameter_dat(data = signaling_responses)

  expect_equal(parameters$rho, parameter_org$rho)
  expect_equal(parameters$bic, parameter_org$bic)
})

test_that("optimize_parameter_hybrid", {

  # read proteomic response file
  signaling_responses <- read.csv(system.file("test_data", "TCGA-BRCA-L4.csv", package = "targetscore"), row.names = 1)

  # Read in biology knowledge base protein interaction
  prior_org <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
    package = "targetscore"
  ))

  # test-data hybrid parameter output
  parameter_org <- readRDS(system.file("test_data_files", "optimize_parameter_hybrid_output.rds",
    package = "targetscore"
  ))
  parameters <- targetscore::optimize_parameter_hybrid(data = signaling_responses, prior = prior_org$wk)

  # Hack for small numeric variations
  rownames(parameters$bic) <- 1:nrow(parameters$bic)
  colnames(parameters$bic) <- 1:ncol(parameters$bic)
  rownames(parameter_org$bic) <- 1:nrow(parameter_org$bic)
  colnames(parameter_org$bic) <- 1:ncol(parameter_org$bic)
  
  # Test equality
  expect_equal(parameters$rho_m, parameter_org$rho_m)
  expect_equal(parameters$rho, parameter_org$rho)
  expect_equal(parameters$kappa, parameter_org$kappa)
  expect_equal(parameters$bic, parameter_org$bic)
})

test_that("get_target_score", {
  # Target Score output
  ts_org <- readRDS(system.file("test_data_files", "get_target_score_output.rds",
    package = "targetscore"
  ))


  # read proteomic responce file
  signaling_responses <- read.csv(system.file("test_data", "TCGA-BRCA-L4.csv",
    package = "targetscore"
  ), row.names = 1)

  # Read in biology knowledge base protein interaction
  network <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
    package = "targetscore"
  ))

  # read proteomic response file
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "targetscore"),
    row.names = 1
  )

  # read functional score value
  fs <- readRDS(system.file("test_data_files", "get_fs_vals_output.rds",
    package = "targetscore"
  ))

  # Calculate Target Score

  ts <- array(0, dim = c(dim(proteomic_responses)[1], dim(proteomic_responses)[2]))
  ts_p <- array(0, dim = c(dim(proteomic_responses)[1], dim(proteomic_responses)[2]))
  ts_q <- array(0, dim = c(dim(proteomic_responses)[1], dim(proteomic_responses)[2]))

  for (i in seq_len(2)) {
    results <- targetscore::get_target_score(
      wk = network$wk,
      wks = network$wks,
      dist_ind = network$dist_ind,
      inter = network$inter,
      n_dose = 1,
      n_prot = dim(proteomic_responses)[2],
      proteomic_responses = proteomic_responses[i, ],
      n_perm = 1,
      verbose = FALSE,
      fs_dat = fs
    )
    ts[i, ] <- results$ts
    ts_p[i, ] <- results$pts
    ts_q[i, ] <- results$q
  }

  colnames(ts) <- colnames(proteomic_responses)
  #ts <- data.frame(rownames(proteomic_responses), ts)
  ts[1,]<-round(ts[1,],2)
  ts[2,]<-round(ts[2,],2)

  colnames(ts_p) <- colnames(proteomic_responses)
  ts_p <- data.frame(rownames(proteomic_responses), ts_p)

  colnames(ts_q) <- colnames(proteomic_responses)
  ts_q <- data.frame(rownames(proteomic_responses), ts_q)

  ts_result <- list(ts = ts, ts_p = ts_p, ts_q = ts_q)
  expect_identical(ts_result$ts, ts_org$ts)
})
