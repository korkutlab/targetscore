#' Get antibody map
#'
#' @return a data.frame with the built in antibody map
#'
#' @concept zeptosens
#' @export
getAntibodyMap <- function() {
  antibodyMapFilename <- system.file("extdata", "antibodyMap.txt", package = "zeptosensUtils")
  antibodyMap <- read.table(antibodyMapFilename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  print(antibodyMap)
  return(antibodyMap)
}
