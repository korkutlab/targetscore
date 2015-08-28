#' Split Sample Names 
#' 
#' @param sampleNames 
#' 
#' @return data.frame a data.frame with the split sample names 
#' 
#' @concept zeptosensPkg
#' @export
splitSampleNames <- function(sampleNames, sampleNameEntries) {
    splitSampleNames <- strsplit(sampleNames, "_")
    
    df <- NULL
    for(i in 1:length(splitSampleNames)) {
        sampleNameNumberIdx <- which(sampleNameEntries == "sampleNameNumber")
        sampleNameIdx <- which(sampleNameEntries == "sampleName")
        treatmentIdx <- which(sampleNameEntries == "treatment")
        timeIdx <- which(sampleNameEntries == "time")
        replicateIdx <- which(sampleNameEntries == "replicate")
        
        # TODO: More generalized 
        tmpDf <- data.frame(
            sampleNameNumber=splitSampleNames[[i]][sampleNameNumberIdx],
            sampleName=splitSampleNames[[i]][sampleNameIdx],
            treatment=splitSampleNames[[i]][treatmentIdx],
            time=splitSampleNames[[i]][timeIdx],
            replicate=splitSampleNames[[i]][replicateIdx],
            stringsAsFactors=FALSE
        )
        
        df <- rbind(df, tmpDf)
    }
    
    return(df)
}