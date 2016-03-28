#' Peak selection
#' 
#' @param arrayData FIXME dataframe from readZeptosensExport
#' @param responseFile FIXME character, name of the output file
#' @param sdctff FIXME standard deviation threshold for selecting responders
#' 
#' @return FIXME
#' 
#' @concept zeptosens
#' @export
selectResponders <- function(arrayData, responseFile, sdctff) {
    arrayData$readout <- log2(abs(arrayData$readout))
    print(arrayData$readout)
    
    dev1 <- 100
    
    colresponders <- paste("rank", "antibody", "readout", "stdev", "treatment", sep = "\t")
    write(colresponders, file = responseFile)
    
    for (iter in 1:100) {
        dev0 <- dev1
        print(dev0)
        dev1 <- sd(arrayData$readout, na.rm = T)
        print(dev1)
        
        if ((dev0 - dev1) < 0.01) {
            break
        }
        responders <- arrayData[which(abs(arrayData$readout) > sdctff * dev1), ]$antibody
        respondrdout <- arrayData[which(abs(arrayData$readout) > sdctff * dev1), 
            ]$readout
        respcond <- arrayData[which(abs(arrayData$readout) > sdctff * dev1), ]$treatment
        arrayData <- arrayData[which(abs(arrayData$readout) < sdctff * dev1), ]
        # print(responders)
        response <- cbind(responders)
        
        for (k in 1:length(responders)) {
            write(paste(iter, response[k], respondrdout[k], dev1, respcond[k], sep = "\t"), 
                file = responseFile, append = T)
        }
        # print(unique(responders))
    }
} 
