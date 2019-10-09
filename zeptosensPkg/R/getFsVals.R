#' Extracted functional score value from COMIC/ONCODRIVE Database. Can be override with Manually set functional score.
#'
#' @param nProt number of proteins tested within the data.
#' @param proteomicResponses Proteomic responses as drug pertubation.
#' @param antibodyMapFile a listing of antibodies, their associated genes, 
#' and modification sites
#' @param fsValueFile a listing of functional scores for each gene manually set up 
#' for overriding COSMIC Database given value, the modification path. (.txt)
#' @param verbose Default as FALSE. If given TRUE, will print out the gene seq mapped with Antibody Map File.
#' @examples
#'
#' @concept zeptosensPkg
#' @export
getFsVals <- function(nProt, proteomicResponses, antibodyMapFile = NULL, fsValueFile = NULL, verbose = F) {

  # match Ab names to gene names & posttranslational modifications
  if (is.null(antibodyMapFile)) {
    antibodyMapFile <- system.file("targetScoreData", "antibodyMap.txt", package = "zeptosensPkg")
  }
  mabToGenes <- antibodyMapFile

  if (verbose) {
    print(mabToGenes)
  }

  if (nProt != ncol(proteomicResponses)) {
    stop("ERROR: nProt is not equal to proteomicResponses column number")
  }

  # Match the protein names in the proteomicresponce with the AntibodyMapfile
  idxAbMap <- which(mabToGenes[, 1] %in% colnames(proteomicResponses))
  if (length(idxAbMap) < nProt) {
    stop("ERROR: Not all columns in data were matched in antibody map")
  }
  if (length(unique(mabToGenes[idxAbMap, 1])) != nProt) {
    print(unique(mabToGenes[idxAbMap, 1]))
    stop("ERROR: Mismatch in the number of selected antibodies and the number of proteomic responses")
  }
  antibodyMapSubset <- mabToGenes[idxAbMap, ]

  mabGenes <- mabToGenes[idxAbMap, 4]
  names(mabGenes) <- mabToGenes[idxAbMap, 1]

  # Get FS value
  mabValue <- mabToGenes[idxAbMap, 6]
  mabFS <- ifelse(mabValue == "a", 1, ifelse(mabValue == "i", -1, ifelse(mabValue == "c", 1, 0)))

  cancerRole <- read.table(system.file("extdata", "Cosmic.txt", package = "zeptosensPkg"), 
                           sep = "\t", header = TRUE, fill = TRUE)
  cosFS <- cancerRole[mabGenes, ]$fs

  fsValue <- mabFS * cosFS

  fs <- data.frame(prot = mabToGenes[idxAbMap, 1], fs = fsValue)

  # Uniqueness of fs value
  prot <- as.character(unlist(unique(fs$prot)))
  fs <- unique(fs)


  fsOverride <- fsValueFile
  # Override with Self setting/external fs value
  if (!is.null(fsOverride)) {
    index <- which(fs$prot %in% fsOverride$prot)
    fs[index, 2] <- fsOverride$fs
  }

  # match with antibody seq
  index <- match(colnames(proteomicResponses), fs$prot)
  fs <- fs[index, ]

  return(fs)
}
