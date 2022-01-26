#' Check the Presence of Antibodies in the TargetScore Antibody Map
#' 
#' @param antibody_vec a vector of antibody labels
#' @param antibody_map if NULL then a default map will be used from antibody_map_file
#' @param antibody_map_file uses an antibody map file from inst/target_score_data 
#' 
#' @details The function will make several edits before the comparison: 
#' * Put an 'X' in front of labels that start with numbers
#' * Replace '-' with '.'
#' * Convert to lower case 
#' @md
#' 
#' @examples 
#' antibody_vec <- readLines(system.file("test_data", "test_antibody_names.txt", package="targetscore"))
#' # BADANTIBODY should be returned as the only entry of 50
#' check_antibody_map_presence(antibody_vec) 
#' 
#' @return a vector of antibodies not present in the TargetScore Antibody Map
#' @export 
check_antibody_map_presence <- function(antibody_vec, 
                                        antibody_map=NULL, 
                                        antibody_map_file = system.file("target_score_data", "antibody_map_12162021.csv", package="targetscore")) {
  # cppa_file <- "~/Downloads/cppa_antibody_list_20211215.txt"
  # tmp <- readLines(cppa_file)
  # antibody_vec <- tmp
  # antibody_map <- read.csv(antibody_map_file)
  
  if(is.null(antibody_map)) {
    antibody_map <- read.csv(antibody_map_file)
  }
  
  antibody_vec_edit <- antibody_vec
  num_start_idx <- grepl("^\\d", antibody_vec_edit)
  antibody_vec_edit[num_start_idx] <- paste0("X", antibody_vec_edit[num_start_idx])
  antibody_vec_edit[num_start_idx]
  
  antibody_vec_edit <- gsub("-", ".", antibody_vec_edit)
  
  idx_present <- which(tolower(antibody_vec_edit) %in% tolower(antibody_map$AntibodyLabel))
  idx_missing <- which(!(tolower(antibody_vec_edit) %in% tolower(antibody_map$AntibodyLabel)))
  
  antibody_vec_edit[idx_missing]  
}
