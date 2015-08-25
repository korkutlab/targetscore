#' Read Zeptosens Export File 
#' 
#' @param inputFile TBA
#' @param sampleNameEntries TBA
#' 
#' @export 
readZeptosensExport <- function(inputFile=NULL, sampleNameEntries=NULL) {
    message("NOTE: Sample Name in *_RFI_Export_Table.xls is assumed to be in the format sampleNumber_sampleName_drug_time_dose_replicateNumber. Use 'samplenNameEntries' parameter to modify as described in documentation.")
    
    colClasses <- c("character", "numeric", "numeric", "character", "character", "character", 
      "numeric", "character", "character", "character", "character", "character", 
      "numeric", "numeric", "numeric", "character", "numeric", "numeric", 
      "numeric", "numeric", "numeric", "numeric", "numeric", "character", 
      "character", "character", "numeric")
    
    data <- read.xlsx2(inputFile, sheetIndex=1, colClasses=colClasses, stringsAsFactors=FALSE)

    splitSampleNames <- strsplit(data[, "Sample.Name"], "_")
    
    df <- NULL
    for(i in 1:nrow(data)) {
        sampleNameNumberIdx <- which(sampleNameEntries == "sampleNameNumber")
        sampleNameIdx <- which(sampleNameEntries == "sampleName")
        treatmentIdx <- which(sampleNameEntries == "treatment")
        timeIdx <- which(sampleNameEntries == "time")
        replicateIdx <- which(sampleNameEntries == "replicate")
        
        tmpDf <- data.frame(
            sampleNumber=data[i, "Sample."],
            sampleNameNumber=splitSampleNames[[i]][sampleNameNumberIdx],
            sampleName=splitSampleNames[[i]][sampleNameIdx],
            treatment=splitSampleNames[[i]][treatmentIdx],
            time=splitSampleNames[[i]][timeIdx],
            replicate=splitSampleNames[[i]][replicateIdx],
            antibody=data[i, "Analyte"],
            readout=data[i, "RFI"],
            cv=data[i, "RFI.CV"],
            quality=data[i, "Class."]
        )
            
        df <- rbind(df, tmpDf)
    }
    
    return(df)
}
