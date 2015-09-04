#' Generate RPPA Graphs
#' 
#' @param platformFile Name of the antibody reference file
#' @param idColumn Column name of IDs
#' @param symbolsColumn Column name of gene symbols
#' @param sitesColumn Column name for phosphorylation sites
#' @param effectColumn Column name for effect of the site on activity
#' @param valuesFile Name of the measurements file
#' @param valueColumn Name of the values column in the measurements file
#' @param valueThreshold The value threshold to be considered as significant
#' @param graphType Either "compatible" or "conflicting"
#' @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
#' 
#' @concept zeptosensPkg
#' @export
generateRPPAGraphs <- function(platformFile, 
                               idColumn="NodeName", 
                               symbolsColumn="Symbol", 
                               sitesColumn="Sites", 
                               effectColumn="Effect", 
                               valuesFile, 
                               valueColumn="Value", 
                               valueThreshold, 
                               graphType="compatible", 
                               outputFilePrefix) {
    
    command <- "generateRPPAGraphs"
    commandJStr <- .jnew("java/lang/String", command)
    
    platformFileJStr <- .jnew("java/lang/String", platformFile)
    idColumnJStr <- .jnew("java/lang/String", idColumn)
    symbolsColumnJStr <- .jnew("java/lang/String", symbolsColumn)
    sitesColumnJStr <- .jnew("java/lang/String", sitesColumn)
    effectColumnJStr <- .jnew("java/lang/String", effectColumn)
    valuesFileJStr <- .jnew("java/lang/String", valuesFile)
    valueColumnJStr <- .jnew("java/lang/String", valueColumn)
    
    valueThresholdJStr <- .jnew("java/lang/Double", valueThreshold)
    graphTypeJStr <- .jnew("java/lang/String", graphType)
    outputFilePrefixJStr <- .jnew("java/lang/String", outputFilePrefix)
    
    argsList <- list(commandJStr, 
                     platformFileJStr, 
                     idColumnJStr, 
                     symbolsColumnJStr, 
                     sitesColumnJStr, 
                     effectColumnJStr,
                     valuesFileJStr,
                     valueColumnJStr,
                     valueThresholdJStr,
                     graphTypeJStr, 
                     outputFilePrefixJStr) 

    .jcall("org/cbio/causality/rppa/RPPAFrontFace","V", command,
           platformFileJStr, 
           idColumnJStr, 
           symbolsColumnJStr, 
           sitesColumnJStr, 
           effectColumnJStr,
           valuesFileJStr,
           valueColumnJStr,
           valueThreshold,
           graphTypeJStr, 
           outputFilePrefixJStr)
    .jcheck()
    
    sifFile <- paste0(outputFilePrefix, ".sif")
    sifCols <- c("PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B", "URIS", "SITES")

    formatFile <- paste0(outputFilePrefix, ".format")
    formatCols <- c("componentType", "componentLabel", "componentProperty", "rgbColor")
    
    sif <- read.table(sifFile, sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE, col.names=sifCols)
    format <- read.table(formatFile, sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE, col.names=formatCols)
    
    results <- list(sif=sif,
                    format=format)
    
    return(results)
}
