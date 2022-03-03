#' Predict Network from signedPC (curated knowledge bio-inferred network)
#'
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#' @param max_dist maximum distance between two antibody. (Default at 1)
#' @param mab_to_genes A list of antibodies, their associated genes, modification sites and effect.
#' @param dist_file A distance file an edgelist with a third column which is the network distance
#'   between the genes in the interaction
#' @param network_file A path to a plain-text, tab-delimited file with the columns: 
#' "PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B", "IDS", "SITES". 
#' "INTERACTION_TYPE" is one of "dephosphorylates", "phosphorylates", "downregulates-expression", or "upregulates-expression"
#' "IDS" is a string of the source of the interaction (e.g., "KEGG")
#' "SITES is a semi-colon separated list of sites that are being phosphorylation or dephosphorylated (e.g., Y398;Y362) otherwise empty string ""
#' @param verbose whether to show debugging information
#' 
#' @note proteomic_responses is used only to retrieve the desired list of 
#' entries for the resulting network
#'
#' @return a list is returned with the following entries:
#' * "wk" inferred network matrix form with edge strength value default at 1 for
#' upregulate and -1 for down regulate.
#' * "wks" as inferred network matrix form with edge strength value default at 1 for
#' upregulate and -1 for down regulate and 2 for phosphorylation and -2 for dephosphorylation.
#' * "dist_ind" A distance file of edgelist with a third column as the network distance
#' between the genes in the interaction.
#' * "inter" file as edgelist of inferred network
#' * "edgelist" file as sif file of edgelist for inferred network.
#'
#' @examples 
#' # READ ANTIBODY FILE ----
#' file <- system.file("target_score_data", "antibody_map_08272020.csv", package = "targetscore")
#' mab_to_genes <- read.csv(file,
#' header = TRUE,
#' stringsAsFactors = FALSE
#' )
#' 
#' # Read proteomic response for cellline1
#' file <- system.file("test_data", "BT474.csv", package = "targetscore")
#' proteomic_responses <- read.csv(file, row.names = 1)
#' 
#'  # Extract network
#'  network <-  predict_bio_network(
#'  n_prot = dim(proteomic_responses)[2],
#'  proteomic_responses = proteomic_responses,
#'  max_dist = 1,
#'  mab_to_genes = mab_to_genes
#'  )
#' 
#' @importFrom utils write.table
#'
#' @concept targetscore
#' @export
predict_bio_network <- function(n_prot, proteomic_responses, max_dist,
                                mab_to_genes, 
                                dist_file = NULL, 
                                network_file=system.file("extdata", "filteredSignedPc_20191113.txt", package = "targetscore"),
                                verbose = FALSE) {
  if (verbose) {
    print(mab_to_genes)
  }

  # pathway distance matrix
  ####################### FIXME##########################
  # if (is.null(dist_file)) {
  #  dist_file <- system.file("target_score_data", "distances.txt", package = "targetscore")
  # }
  # tmp_dist <- read.table(dist_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  ######################################################

  if (n_prot != ncol(proteomic_responses)) {
    stop("ERROR: n_prot is not equal to proteomic_responses column number")
  }

  # Filter dist to only keep those with a distance less than max_dist
  # idx <- which(tmp_dist[, 3] <= max_dist)
  # dist <- tmp_dist[idx, ]
  # dist
  idx_ab_map <- which(mab_to_genes[, 1] %in% colnames(proteomic_responses))
  #    print(mab_to_genes[, 1])
  #    print(colnames(proteomic_responses))
  #    print(unique(mab_to_genes[idx_ab_map,]))
  #    if(length(idx_ab_map) < n_prot) {
  #    print(length(unique(mab_to_genes[idx_ab_map,1])))
  if (length(idx_ab_map) < n_prot) {
    #        print(length(idx_ab_map))
    stop("ERROR: Not all columns in data were matched in antibody map")
  }
  #    print((unique(mab_to_genes[idx_ab_map, 1])))
  if (length(unique(mab_to_genes[idx_ab_map, 1])) != n_prot) {
    print(unique(mab_to_genes[idx_ab_map, 1]))
    stop("ERROR: Mismatch in the number of selected antibodies and the number of proteomic responses")
  }

  # Used by match_genes_to_edgelist to account for the cases where multiple entries for the same antibody
  # exist in the antibody map
  # antibody_map_subset <- mab_to_genes[idx_ab_map, ]

  mab_genes <- mab_to_genes[idx_ab_map, 4]
  names(mab_genes) <- mab_to_genes[idx_ab_map, 1]
  ########################### FIXME#################################################
  # dist_list <- match_genes_to_edgelist(
  #  genes1 = mab_genes, genes2 = NULL, annot_edgelist = dist,
  #  antibody_vec = colnames(proteomic_responses), use_annot = TRUE, verbose = TRUE
  # )
  #################################################################################
  # # Get interactions for measured genes
  # dist_gene1 <- pmatch(dist[, 1], mab_to_genes[measured_genes, 4], duplicates.ok = TRUE)
  # dist_gene1Name <- mab_to_genes[measured_genes, 4][dist_gene1]
  # dist_gene2 <- pmatch(dist[, 2], mab_to_genes[measured_genes, 4], duplicates.ok = TRUE)
  # dist_gene2Name <- mab_to_genes[measured_genes, 4][dist_gene2]
  #
  # # distance framework
  #################### FIXME############################
  # dist_list <- data.frame(dist_gene1=dist_gene1, dist_gene2=dist_gene2, dist=dist[, 3],
  #                         dist_gene1Name=dist_gene1Name, dist_gene2Name=dist_gene2Name, stringsAsFactors = FALSE)
  #####################################################
  # dist_list[is.na(dist_list[, 1]), 3] <- 100
  # dist_list[is.na(dist_list[, 2]), 3] <- 100
  # dist_list[, 1]

  # dist_ind(upstream,downstream)
  dist_ind <- matrix(Inf,
    ncol = n_prot, nrow = n_prot,
    dimnames = list(colnames(proteomic_responses), colnames(proteomic_responses))
  )
  #### Distance File not functional #########FIXME
  # for (i in 1:length(dist_list[, 1])) {
  #  dist_ind[dist_list[i, 1], dist_list[i, 2]] <- dist_list[i, 3]

  #    if (dist_ind[dist_list[i, 1], dist_list[i, 2]] > max_dist) {
  #    dist_ind[dist_list[i, 1], dist_list[i, 2]] <- Inf
  #  }

  #    if (dist_ind[dist_list[i, 1], dist_list[i, 2]] == 0) {
  #    dist_ind[dist_list[i, 1], dist_list[i, 2]] <- Inf
  # }
  # }
  ##################################################

  ### get the network product### phospFile <- system.file('SignedPC', 'phosphorylates.txt',
  ### package='targetscore') phosp <- read.csv(phospFile, sep='\t', header=TRUE, na.strings =
  ### c('', ' ')) phosp <- phosp[,-3] dephospFile <- system.file('SignedPC',
  ### 'dephosphorylates.txt', package='targetscore') dephosp <- read.csv(dephospFile, sep='\t',
  ### header=TRUE, na.strings = c('', ' ')) dephosp <- dephosp[,-3] upexpFile <-
  ### system.file('SignedPC', 'dephosphorylates.txt', package='targetscore') upexp <-
  ### read.csv(upexpFile, sep='\t', header=TRUE, na.strings = c('', ' ')) upexp <- upexp[,-3]
  ### dwnexpFile <- system.file('SignedPC', 'downregulates-expression.txt', package='targetscore')
  ### dwnexp <- read.csv(dwnexpFile, sep='\t', header=TRUE, na.strings = c('', ' ')) dwnexp <-
  ### dwnexp[,-3]

  #    results <- downloadSignedPC(forceCache=TRUE)
  #    results <- read.table(file="results_network1.txt",header=T)
  results <- read.table(network_file, sep = "\t", header = TRUE, 
                        fill = TRUE, stringsAsFactors = FALSE)
  # write.table(results, file="results_network1.txt",quote=F)
  
  dephosp <- results[which(results$INTERACTION_TYPE == "dephosphorylates"),]
  phosp <- results[which(results$INTERACTION_TYPE == "phosphorylates"),]
  dwnexp <- results[which(results$INTERACTION_TYPE == "downregulates-expression"),]
  upexp <- results[which(results$INTERACTION_TYPE == "upregulates-expression"),]

  # NOTE: SIF has interaction type as column 2, edgelists (like distances) do
  # not have this, so convert the SIF to an edgelist
  dephosp <- dephosp[, c(1, 3)]
  phosp <- phosp[, c(1, 3)]
  dwnexp <- dwnexp[, c(1, 3)]
  upexp <- upexp[, c(1, 3)]

  # only concentration nodes are included in up & downregulation
  # mabToGenes_c <- mab_to_genes[which(mab_to_genes$Effect == "c"), ]

  # define wk
  wk <- matrix(0,
    ncol = n_prot, nrow = n_prot,
    dimnames = list(colnames(proteomic_responses), colnames(proteomic_responses))
  )
  wks <- matrix(0,
    ncol = n_prot, nrow = n_prot,
    dimnames = list(colnames(proteomic_responses), colnames(proteomic_responses))
  )

  # upregulation expression, wk=1
  # upexp_gene1 <- pmatch(upexp[, 1], mabToGenes_c[measured_genes, 4], duplicates.ok = TRUE)
  # upexp_gene2 <- pmatch(upexp[, 3], mabToGenes_c[measured_genes, 4], duplicates.ok = TRUE)
  # upexp_gene <- cbind(upexp_gene1, upexp_gene2)

  # Define genes by their effects
  tmp_idx_c <- intersect(idx_ab_map, which(mab_to_genes$Effect == "c"))
  tmp_idx_ac <- intersect(idx_ab_map, which(mab_to_genes$Effect != "i"))
  tmp_idx_ai <- intersect(idx_ab_map, which(mab_to_genes$Effect != "c"))
  tmp_idx_a <- intersect(idx_ab_map, which(mab_to_genes$Effect == "a"))

  tmp_genes_c <- mab_to_genes[tmp_idx_c, 4]
  tmp_genes_a <- mab_to_genes[tmp_idx_ac, 4]
  tmp_genes_d <- mab_to_genes[tmp_idx_ai, 4]
  tmp_genes_ao <- mab_to_genes[tmp_idx_a, 4]

  names(tmp_genes_c) <- mab_to_genes[tmp_idx_c, 1]
  names(tmp_genes_a) <- mab_to_genes[tmp_idx_ac, 1]
  names(tmp_genes_d) <- mab_to_genes[tmp_idx_ai, 1]
  names(tmp_genes_ao) <- mab_to_genes[tmp_idx_a, 1]

  # saveRDS(colnames(proteomic_responses), "colnames_proteomic_responses.rds")

  # only concentration and act. phospho nodes are included in up & downregulation
  upexp_gene <-  match_genes_to_edgelist(
    genes1 = tmp_genes_a, genes2 = tmp_genes_c, annot_edgelist = upexp,
    antibody_vec = colnames(proteomic_responses), use_annot = FALSE, verbose = verbose
  )

  # saveRDS(upexp_gene, "upexp_gene.rds")

  for (i in 1:length(upexp[, 1])) {
    wk[upexp_gene[i, 1], upexp_gene[i, 2]] <- 1
    wks[upexp_gene[i, 1], upexp_gene[i, 2]] <- 1

    if (verbose) {
      cat(upexp_gene[i, 1], "\n")
    }
  }

  # downregulation expression, wk=-1
  dwnexp_gene <-  match_genes_to_edgelist(
    genes1 = tmp_genes_a, genes2 = tmp_genes_c, annot_edgelist = dwnexp,
    antibody_vec = colnames(proteomic_responses), use_annot = FALSE, verbose = verbose
  )
  # cov318 results in 15

  for (i in 1:length(dwnexp[, 1])) {
    wk[dwnexp_gene[i, 1], dwnexp_gene[i, 2]] <- -1
    wks[dwnexp_gene[i, 1], dwnexp_gene[i, 2]] <- -1
  }

  # phosphorylates wk=1 only active and concentration states are upstream
  # mabToGenes_a <- mab_to_genes[which(mab_to_genes$Effect != "i"), ]
  # mabToGenes_d <- mab_to_genes[which(mab_to_genes$Sites != "c"), ]
  # phos_gene1 <- pmatch(phosp[, 1], mabToGenes_a[measured_genes, 4], duplicates.ok = TRUE)
  # phos_gene2 <- pmatch(phosp[, 3], mabToGenes_d[measured_genes, 4], duplicates.ok = TRUE)
  # phos_gene <- cbind(phos_gene1, phos_gene2)

  phos_gene <-  match_genes_to_edgelist(
    genes1 = tmp_genes_ao, genes2 = tmp_genes_d, annot_edgelist = phosp,
    antibody_vec = colnames(proteomic_responses), use_annot = FALSE, verbose = verbose
  )
  # cov318 13 results

  for (i in 1:length(phos_gene[, 1])) {
    wk[phos_gene[i, 1], phos_gene[i, 2]] <- 1
    wks[phos_gene[i, 1], phos_gene[i, 2]] <- 2
  }

  # dephosphorylates wk=-1 only active and concentration states are upstream
  # mabToGenes_a <- mab_to_genes[which(mab_to_genes$Effect != "i"), ]
  # mabToGenes_d <- mab_to_genes[which(mab_to_genes$Sites != "c"), ]
  # dephos_gene1 <- pmatch(dephosp[, 1], mabToGenes_a[measured_genes, 4], duplicates.ok = TRUE)
  # dephos_gene2 <- pmatch(dephosp[, 3], mabToGenes_d[measured_genes, 4], duplicates.ok = TRUE)
  # dephos_gene <- cbind(dephos_gene1, dephos_gene2)

  dephos_gene <-  match_genes_to_edgelist(
    genes1 = tmp_genes_a, genes2 = tmp_genes_d, annot_edgelist = dephosp,
    antibody_vec = colnames(proteomic_responses), use_annot = FALSE, verbose = verbose
  )

  for (i in 1:length(dephos_gene[, 1])) {
    wk[dephos_gene[i, 1], dephos_gene[i, 2]] <- -1
    wks[dephos_gene[i, 1], dephos_gene[i, 2]] <- -2
  }

  inter <- (which(wk != 0, arr.ind = TRUE))

  if (verbose) {
    print(inter)
  }

  # dist_file

  for (i in 1:n_prot) {
    for (j in 1:n_prot) {
      if (wk[i, j] != 0) {
        dist_ind[i, j] <- 1
      } else {
        dist_ind[i, j] <- Inf
      }
    }
  }
  # Network to edgelist
  edgelist <-  create_sif_from_matrix(
    t_net = wk,
    col_genelist = colnames(wk),
    row_genelist = rownames(wk)
  )


  network <- list(wk = wk, wks = wks, dist_ind = dist_ind, inter = inter, edgelist = edgelist)

  return(network)
}
