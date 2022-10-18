#' Extracted functional score value from COMIC database; this can be overriden with manually curated scores
#'
#' @param n_prot Antibody number of input data.
#' @param proteomic_responses Input drug perturbation data. With columns as antibody, rows as samples.
#' @param mab_to_genes A list of antibodies, their associated genes, modification sites and effect.
#' @param fs_override a listing of functional scores for each gene manually set up
#' for overriding COSMIC Database given value, the modification path. (.txt)
#' @param cancer_role_file File with at two least columns: "Gene.Symbol" and 
#'   "Role.in.Cancer" where the "Role.in.Cancer" has mainly values "TSG" 
#'   (tumor suppressor gene) or "oncogene"; "TSG, fusion" and "oncogene, fusion" 
#'   are also understood. Other values (e.g., "TSG, oncogene", "NA", or blank) 
#'   will be interpreted as dual or unknown as part of the functional score. 
#'   One source for this is the COSMIC 
#'   Cancer Gene Census: https://cancer.sanger.ac.uk/cosmic/download
#' @param verbose Default as FALSE. If given TRUE, will print out the gene seq mapped with antibody map file.
#' 
#' @return 
#' * "fs_final" dataframe with two coloumns: prot as antibody label; fs as functional #' score 
#'
#' @examples
#' # Read fs_manually set file
#' fs_override_org <- readRDS(system.file("test_data_files", "fs_value_file.rds",
#' package = "targetscore"
#' ))
#' 
#' # Read proteomic response file
#' file <- system.file("test_data", "BT474.csv", package = "targetscore")
#' proteomic_responses <- read.csv(file, row.names = 1)
#' 
#' # Read antibody file
#' file <- system.file("target_score_data", "antibody_map.csv", package = "targetscore")
#' mab_to_genes <- read.csv(file,
#' header = TRUE,
#' stringsAsFactors = FALSE
#' )
#' 
#' fs <- get_fs_vals(
#' n_prot = ncol(proteomic_responses), proteomic_responses = proteomic_responses,
#' mab_to_genes = mab_to_genes, fs_override = fs_override_org
#' )
#' 
#' @importFrom utils read.table
#'
#' @concept targetscore
#' @export
get_fs_vals <- function(n_prot, proteomic_responses, mab_to_genes, fs_override=NULL, 
  cancer_role_file=system.file("extdata", "cancer_gene_census_v88.csv", package="targetscore"), verbose=FALSE) {
  
  if (verbose) {
    print(mab_to_genes)
  }

  if (n_prot != ncol(proteomic_responses)) {
    stop("ERROR: n_prot is not equal to proteomic_responses column number")
  }

  # Match the protein names in the proteomic responce with the antibody map file
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

  # Get Functional Score Value
  ## Get the protein/phosphoprotein effect information from the antibody map
  mab_value <- mab_to_genes[idx_ab_map, 6]
  mab_phospho_effect <- ifelse(mab_value == "a", 1, ifelse(mab_value == "i", -1, ifelse(mab_value == "c", 1, 0)))
  
  cancer_role <- read.csv(cancer_role_file)
  cancer_role <- cancer_role[, c("Gene.Symbol", "Role.in.Cancer")]
  
  # For consistency with existing code
  colnames(cancer_role) <- c("gene", "role_in_cancer")
  
  # Mapping: oncogene as +1, tumor suppressor as -1, or 0 both or unknown
  cancer_role$fs <- 0
  cancer_role$fs[cancer_role$role_in_cancer == "oncogene" | cancer_role$role_in_cancer == "oncogene, fusion"] <- 1
  cancer_role$fs[cancer_role$role_in_cancer == "TSG" | cancer_role$role_in_cancer == "TSG, fusion"] <- -1
  
  index <- match(mab_genes, cancer_role$gene)
  cos_fs <- cancer_role[index, ]$fs
  cos_fs[is.na(cos_fs)] <- 0

  ## Final functional score values will the multiplication of the protein effect and the functional annotation
  fs_value <- mab_phospho_effect * cos_fs

  fs <- data.frame(prot = mab_to_genes[idx_ab_map, 1], fs = fs_value, stringsAsFactors = FALSE)
  
  # Replace NAs with 0 (i.e., unknown)
  fs$fs[is.na(fs$fs)] <- 0
  
  # Uniqueness of fs value
  fs <- unique(fs)

  # Deal with duplicate and double value extracted
  prot_dup <- fs$prot[duplicated(fs$prot)]
  
  if(length(prot_dup) > 0) {
    for (i in 1:length(prot_dup)) {
      if(verbose) {
        message("DEBUG: Duplicate: ",  prot_dup[i], "\n")
      }
      
      index <- which(fs$prot %in% prot_dup[i])
      fs_dup <- fs[index, ]
      
      # Remove the duplicate
      fs <- fs[-index, ]
      
      # Replace the former duplicates with a new entry based on these rules
      if (any(c(fs_dup$fs) == 0) & any(c(fs_dup$fs) == 1)) {
        dup <- c(as.character(prot_dup[i]), 1)
      }
      if (any(c(fs_dup$fs) == 1) & any(c(fs_dup$fs) == -1)) {
        dup <- c(as.character(prot_dup[i]), 0)
      }
      if (any(c(fs_dup$fs) == 0) & any(c(fs_dup$fs) == -1)) {
        dup <- c(as.character(prot_dup[i]), -1)
      }
      if (any(c(fs_dup$fs) == 0) & any(c(fs_dup$fs) == -1) & any(c(fs_dup$fs) == 1)) {
        dup <- c(as.character(prot_dup[i]), 0)
      }
      
      fs <- rbind(fs, dup)
    }    
  }

  # Override with self setting/external fs value
  if (!is.null(fs_override)) {
    # Select fs_override
    index <- match(fs$prot, fs_override$prot)
    fs_override <- fs_override[index, ]

    # override fs with fs_override
    index <- match(fs_override$prot, fs$prot)
    fs[index, 2] <- fs_override$fs
  }

  # match with antibody seq
  index <- match(colnames(proteomic_responses), fs$prot)
  fs <- fs[index, ]

  # as.numeric
  fs_final <- data.frame(prot = fs$prot, fs = as.numeric(fs$fs), stringsAsFactors = FALSE)
  fs_final <- fs_final[order(fs_final$prot),]
  
  return(fs_final)
}
