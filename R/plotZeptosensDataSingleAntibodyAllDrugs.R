#' Plot Zeptosens Data
#' 
#' @param normArray1 normalized data from array 1
#' @param normArray2 normalized data from array 2
#' @param antibodies vector of antibodies 
#' @param excludeTreatments vector of treatments to be excluded
#' @param ctrlTreatment string, treatment label used as control 
#' @param ctrlTime string, time label used as control 
#' @param plotColNames vector, array names 
#' @param plotRowNames vector, treatment time names
#' @param plotColors vector, colors for the two arrays
#' @param plotDir string, directory to put plots
#' 
#' @concept zeptosensPkg
#' @export
plotZeptosensDataSingleAntibodyAllDrugs <- function(normArray1, normArray2, antibodies, excludeTreatments, 
                                                    ctrlTreatment, ctrlTime, plotColNames, plotRowNames, plotColors, plotDir) {
    allDrugs <- unique(normArray1[, "treatment"])
    tmpDrugs <- allDrugs[!(unique(normArray1[, "treatment"]) %in% excludeTreatments)]

    for(curAntibody in antibodies) {
#         plotRows <- ceiling(sqrt(length(antibodies)))
#         plotCols <- ceiling(sqrt(length(antibodies)))
        
        tmpAntibody <- gsub("[[:punct:]]", "-", curAntibody)
        pdf(paste0(plotDir, tmpAntibody, ".pdf"), width=20, height=8.5)
        
#        par(mfrow=c(plotRows, plotCols))
        
        #par(mar=c(2,2,1.5,1), oma=c(0,0,2,0))
        
        # Counter for colors
        cnt <- 1
        
        # Setup plot
        plot(c(1, length(tmpDrugs)), range(normArray1$readout), main=paste0(curAntibody, " Sample: ", plotColNames), 
             xaxt="n", xlab="Time", ylab="RFI") 
        
        axis(1, at=1:length(tmpDrugs), labels=plotRowNames)
        
        tmpMat <- NULL
        
        for(curDrug in tmpDrugs) {
            tmpArray1Entries <- normArray1[which(normArray1[, "antibody"] == curAntibody), ]
#            tmpArray2Entries <- normArray2[which(normArray2[, "antibody"] == curAntibody), ]
            
            idx <- intersect(which(tmpArray1Entries[, "antibody"] == curAntibody), which(tmpArray1Entries[, "treatment"] == curDrug))
            
            t0Idx <- intersect(which(tmpArray1Entries[, "treatment"] == ctrlTreatment), which(tmpArray1Entries[, "time"] == ctrlTime))
            allIdx <- c(t0Idx, idx)

            tmpMat <- cbind(tmpMat, tmpArray1Entries[allIdx, "normReadout"])
            
            lines(1:length(allIdx), tmpArray1Entries[allIdx, "normReadout"], type="b", lwd=1.5, col=plotColors[cnt]) 
            
            cnt <- cnt + 1
        }
        
        colnames(tmpMat) <- tmpDrugs
        rownames(tmpMat) <- plotRowNames
        
        tmpMat <- t(tmpMat)
        
        #plot(tmpMat, main=curAntibody, col="blue")                
        #title(main=curAntibody, outer=TRUE)
        
        legend(2, 2, # Position
               tmpDrugs, # Text
               lty=rep(1, length(tmpDrugs)), # Legend symbols
               lwd=rep(2.4, length(tmpDrugs)), # Symbol width
               col=plotColors) # Colors
        
        dev.off()
    }
    
    return(tmpMat)
}