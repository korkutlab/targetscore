#' Generate Causality Graph
#' 
#' @param platformFile Name of the antibody reference file
#' @param idColumn Column name of IDs
#' @param symbolsColumn Column name of gene symbols
#' @param sitesColumn Column name for phosphorylation sites
#' @param effectColumn Column name for effect of the site on activity
#' @param valuesFile Name of the measurements file
#' @param valueColumn Name of the values column in the measurements file
#' @param valueThreshold The value threshold to be considered as significant
#' @param graphType Either 'compatible' or 'conflicting'
#' @param doSiteMatch option to enforce matching a phosphorylation site in the network with the annotation of antibody
#' @param siteMatchProximityThreshold when site matching is on, this parameter sets the proxomity threshold for a
#'   site number in the relation to match the site whose change is observed in the data (DEFAULT: 0)
#' @param siteEffectProximityThreshold when the site effect is not known, we can approximate it with the known
#'   effect of proximate sites. This parameter sets the proximity threshold for using the proximate sites for that 
#'   prediction. (DEFAULT: 0)
#' @param geneCentric A boolean option to produce a gene-centric (TRUE) or an antibody-centric (FALSE) graph
#' @param colorSaturationValue The value that maps to the most saturated color (DEFAULT: 10)
#' @param outputFilePrefix If the user provides xxx, then xxx.sif and xxx.format are generated
#' @param customNetworkDirectory The directory that the network will be downloaded and SignedPC
#'   directory will be created in. Pass null to use default.
#' 
#' @concept zeptosensPkg
#' @export
generateCausalityGraph <- function(platformFile, idColumn = "NodeName", symbolsColumn = "Symbol", sitesColumn = "Sites", 
    effectColumn = "Effect", valuesFile, valueColumn = "Value", valueThreshold, graphType = "compatible", 
    doSiteMatch = TRUE, siteMatchProximityThreshold = 0, siteEffectProximityThreshold = 0, geneCentric = TRUE, 
    colorSaturationValue=10, outputFilePrefix, customNetworkDirectory = tempdir()) {
    
    command <- "generateCausalityGraph"
    commandJStr <- .jnew("java/lang/String", command)
    
    platformFileJStr <- .jnew("java/lang/String", platformFile)
    idColumnJStr <- .jnew("java/lang/String", idColumn)
    symbolsColumnJStr <- .jnew("java/lang/String", symbolsColumn)
    sitesColumnJStr <- .jnew("java/lang/String", sitesColumn)
    effectColumnJStr <- .jnew("java/lang/String", effectColumn)
    valuesFileJStr <- .jnew("java/lang/String", valuesFile)
    valueColumnJStr <- .jnew("java/lang/String", valueColumn)
    graphTypeJStr <- .jnew("java/lang/String", graphType)
    outputFilePrefixJStr <- .jnew("java/lang/String", outputFilePrefix)
    customNetworkDirectoryJStr <- .jnew("java/lang/String", customNetworkDirectory)
    colorSaturationValue <- .jfloat(colorSaturationValue)
    
    .jcall("org/panda/causalpath/run/CausalityAnalysisSingleMethodInterface", "V", command, 
           platformFileJStr, idColumnJStr, symbolsColumnJStr, sitesColumnJStr, effectColumnJStr, 
           valuesFileJStr, valueColumnJStr, valueThreshold, graphTypeJStr, doSiteMatch, 
           as.integer(siteMatchProximityThreshold), as.integer(siteEffectProximityThreshold), 
           geneCentric, colorSaturationValue, outputFilePrefixJStr, customNetworkDirectoryJStr)
    .jcheck()
    
    sifFile <- paste0(outputFilePrefix, ".sif")
    sifCols <- c("PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B", "URIS", "SITES")
    
    formatFile <- paste0(outputFilePrefix, ".format")
    formatCols <- c("componentType", "componentLabel", "componentProperty", "rgbColor")
    
    sif <- read.table(sifFile, sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE, 
        col.names = sifCols)
    format <- read.table(formatFile, sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE, 
        col.names = formatCols)
    
    results <- list(sif = sif, format = format)
    
    return(results)
} 
