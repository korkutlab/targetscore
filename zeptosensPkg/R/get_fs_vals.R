#' Extracted functional score value from COMIC/ONCODRIVE Database. Can be override with Manually set functional score.
#'
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#' @param mab_to_genes A list of antibodies, their associated genes, modification sites and effect.
#' @param fs_value_file a listing of functional scores for each gene manually set up
#' for overriding COSMIC Database given value, the modification path. (.txt)
#' @param verbose Default as FALSE. If given TRUE, will print out the gene seq mapped with Antibody Map File.
#'
#' @importFrom utils read.table
#'
#' @concept zeptosensPkg
#' @export
get_fs_vals <- function(n_prot, proteomic_responses, mab_to_genes, fs_value_file = NULL, verbose = FALSE) {
  if (verbose) {
    print(mab_to_genes)
  }

  if (n_prot != ncol(proteomic_responses)) {
    stop("ERROR: n_prot is not equal to proteomic_responses column number")
  }

  # Match the protein names in the proteomicresponce with the AntibodyMapfile
  idx_ab_map <- which(mab_to_genes[, 1] %in% colnames(proteomic_responses))
  if (length(idx_ab_map) < n_prot) {
    stop("ERROR: Not all columns in data were matched in antibody map")
  }
  if (length(unique(mab_to_genes[idx_ab_map, 1])) != n_prot) {
    print(unique(mab_to_genes[idx_ab_map, 1]))
    stop("ERROR: Mismatch in the number of selected antibodies and the number of proteomic responses")
  }

  mab_genes <- mab_to_genes[idx_ab_map, 4]
  names(mab_genes) <- mab_to_genes[idx_ab_map, 1]

  # Get FS value
  mab_value <- mab_to_genes[idx_ab_map, 6]
  mab_fs <- ifelse(mab_value == "a", 1, ifelse(mab_value == "i", -1, ifelse(mab_value == "c", 1, 0)))

  cancer_role <- read.table(system.file("extdata", "Cosmic.txt", package = "zeptosensPkg"),
    sep = "\t", header = TRUE, fill = TRUE
  )
  cos_fs <- cancer_role[mab_genes, ]$fs

  fs_value <- mab_fs * cos_fs

  fs <- data.frame(prot = mab_to_genes[idx_ab_map, 1], fs = fs_value)

  # Uniqueness of fs value
  fs <- unique(fs)

  fs_override <- fs_value_file
  # Override with Self setting/external fs value
  if (!is.null(fs_override)) {
    index <- which(fs$prot %in% fs_override$prot)
    fs[index, 2] <- fs_override$fs
  }

  # match with antibody seq
  index <- match(colnames(proteomic_responses), fs$prot)
  fs <- fs[index, ]

  return(fs)
}
