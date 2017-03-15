#' Split Sample Names 
#' 
#' @param sampleNames 
#' 
#' @return data.frame a data.frame with the split sample names 
#' 
#' @concept zeptosens
#' @export
splitSampleNames <- function(sampleNames, sampleNameEntries) {
    splitSampleNames <- strsplit(sampleNames, "_")
    
    df <- NULL
    
    for (i in 1:length(splitSampleNames)) {
        sampleNumberIdx <- which(sampleNameEntries == "sNumber")
        sampleIdx <- which(sampleNameEntries == "sample")
        treatmentIdx <- which(sampleNameEntries == "treatment")
        doseIdx <- which(sampleNameEntries == "dose")
        timeIdx <- which(sampleNameEntries == "time")
        replicateIdx <- which(sampleNameEntries == "replicate")
        dateIdx <- which(sampleNameEntries == "date")
        notesIdx <- which(sampleNameEntries == "notes")
        
        availableCols <- c(sNumber = sampleNumberIdx, sample = sampleIdx, treatment = treatmentIdx, 
            dose = doseIdx, time = timeIdx, replicate = replicateIdx, date = dateIdx, 
            notes = notesIdx)
        
        tmpRow <- NULL
        
        for (col in names(availableCols)) {
            if (col == "sNumber") {
                tmpRow <- c(tmpRow, splitSampleNames[[i]][sampleNumberIdx])
            }
            
            if (col == "sample") {
                tmpRow <- c(tmpRow, splitSampleNames[[i]][sampleIdx])
            }
            
            if (col == "treatment") {
                tmpRow <- c(tmpRow, splitSampleNames[[i]][treatmentIdx])
            }
            
            if (col == "dose") {
                tmpRow <- c(tmpRow, splitSampleNames[[i]][doseIdx])
            }
            
            if (col == "time") {
                tmpRow <- c(tmpRow, splitSampleNames[[i]][timeIdx])
            }
            
            if (col == "replicate") {
                tmpRow <- c(tmpRow, splitSampleNames[[i]][replicateIdx])
            }
          
            if (col == "date") {
              tmpRow <- c(tmpRow, splitSampleNames[[i]][dateIdx])
            }
            
            if (col == "notes") {
                tmpRow <- c(tmpRow, splitSampleNames[[i]][notesIdx])
            }
        }
        
        df <- rbind(df, tmpRow)
    }
    
    df <- as.data.frame(df, stringsAsFactors = FALSE, row.names = 1:length(sampleNames))
    colnames(df) <- sampleNameEntries[which(!is.na(sampleNameEntries))]
    
    return(df)
} 
