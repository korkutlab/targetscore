test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("match_genes_to_edgelist", {
  antibody_map_file <- system.file("targetScoreData", "antibodyMap.txt", package = "zeptosensPkg")
  mab_to_genes <- read.table(antibody_map_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  proteomic_responses_file <- system.file("test_data", "BT474.csv", package = "zeptosensPkg")
  proteomic_responses <- read.csv(proteomic_responses_file, row.names = 1)

  dist_file <- system.file("targetScoreData", "distances.txt", package = "zeptosensPkg")
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

  snapshot_file <- system.file("test_data_files", "match_genes_to_edgelist_test_output.rds", package = "zeptosensPkg")
  snapshot <- readRDS(snapshot_file)

  expect_identical(dist_list, snapshot)
})

test_that("predict_bio_network", {
  skip_on_cran()

  # READ ANTIBODY FILE ----
  mab_to_genes <- read.table(system.file("targetscoreData", "antibodyMapFile.txt", package = "zeptosensPkg"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )

  # Read proteomic response for cellline1
  proteomic_responses <- read.csv(system.file("test_data", "BT474.csv", package = "zeptosensPkg"), row.names = 1)

  # Extract network
  network <- zeptosensPkg::predict_bio_network(
    n_prot = dim(proteomic_responses)[2],
    proteomic_responses = proteomic_responses,
    max_dist = 1,
    mab_to_genes = mab_to_genes
  )

  network_org <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
    package = "zeptosensPkg"
  ))

  expect_identical(network$wk, network_org$wk)
  expect_identical(network$wks, network_org$wks)
  expect_identical(network$dist_ind, network_org$dist_ind)
  expect_identical(network$inter, network_org$inter)
})

test_that("create_sif_from_matrix", {
  network_org <- readRDS(system.file("test_data_files", "predict_bio_network_network_output.rds",
    package = "zeptosensPkg"
  ))

  wk_org <- readRDS(system.file("test_data_files", "create_sif_from_matrix_output.rds",
    package = "zeptosensPkg"
  ))

  edgelist_wk <- create_sif_from_matrix(t_net = network_org$wk, genelist = colnames(network_org$wk))

  expect_identical(edgelist_wk, wk_org)
})

test_that("get_fs_vals", {
  network_org <- readRDS(system.file("test_data_files", "get_fs_vals_output.rds",
    package = "zeptosensPkg"
  ))

  wk_org <- readRDS(system.file("test_data_files", "create_sif_from_matrix_output.rds",
    package = "zeptosensPkg"
  ))

  fs <- get_fs_vals()

  expect_identical(edgelist_wk, wk_org)
})
