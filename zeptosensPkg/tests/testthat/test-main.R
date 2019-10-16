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

  expect_equal(dist_list, snapshot)
})
