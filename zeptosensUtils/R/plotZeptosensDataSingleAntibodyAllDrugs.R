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
plotZeptosensDataSingleAntibodyAllDrugs <- function(normArray1, normArray2, 
                                                    antibodies, excludeTreatments, 
                                                    ctrlTreatment, ctrlTime, 
                                                    plotColNames, plotRowNames, 
                                                    plotColors, plotDir, 
                                                    asPdf=TRUE, plotPrefix=NULL) {
    allDrugs <- unique(normArray1[, "treatment"])
    tmpDrugs <- allDrugs[!(unique(normArray1[, "treatment"]) %in% excludeTreatments)]

    for(curAntibody in antibodies) {
        tmpMat <- NULL
        
        for(curDrug in tmpDrugs) {
            tmpArray1Entries <- normArray1[which(normArray1[, "antibody"] == curAntibody), ]
#            tmpArray2Entries <- normArray2[which(normArray2[, "antibody"] == curAntibody), ]
            
            idx <- intersect(which(tmpArray1Entries[, "antibody"] == curAntibody), 
                             which(tmpArray1Entries[, "treatment"] == curDrug))
            
            t0Idx <- intersect(which(tmpArray1Entries[, "treatment"] == ctrlTreatment), 
                               which(tmpArray1Entries[, "time"] == ctrlTime))
            allIdx <- c(t0Idx, idx)

            tmpMat <- cbind(tmpMat, tmpArray1Entries[allIdx, "normReadout"])
        }
        
        colnames(tmpMat) <- tmpDrugs
        rownames(tmpMat) <- plotRowNames

        tmpAntibody <- gsub("[[:punct:]]", "-", curAntibody)
        
        if(!is.null(plotPrefix)) {
            fileName <- paste0(plotDir, plotPrefix, "_", tmpAntibody, ".pdf")
        } else {
            fileName <- paste0(plotDir, tmpAntibody, ".pdf")
        }
        
        pdf(fileName, width=9, height=7)
        
        # par(mfrow=c(plotRows, plotCols))
        # par(mar=c(2,2,1.5,1), oma=c(0,0,2,0))
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        
        # Setup plot
        # type="n" for nothing 
        plot(c(1, length(tmpDrugs)), range(tmpMat), main=paste0(curAntibody, " Sample: ", plotColNames), 
             xaxt="n", xlab="Time", ylab="RFI", type="n") 
        plotChar <- seq(18, 18+length(tmpDrugs), 1)
        
        axis(1, at=1:length(tmpDrugs), labels=plotRowNames)
        
        for(i in 1:length(tmpDrugs)) {
            curDrug <- tmpDrugs[i]
            lines(1:nrow(tmpMat), tmpMat[, curDrug], type="b", lwd=1.5, col=plotColors[i], pch=plotChar[i]) 
        }
        
        #plot(tmpMat, main=curAntibody, col="blue")                
        #title(main=curAntibody, outer=TRUE)
        legend("topright", # Position
               legend=tmpDrugs, # Text
               lty=rep(1, length(tmpDrugs)), # Legend line
               pch=plotChar, # Legend symbol
               lwd=rep(2.4, length(tmpDrugs)), # Symbol width
               col=plotColors, # Colors
               inset=c(-0.2, 0),
               title="Treatment") 
        
        dev.off()
    }
    
    return(tmpMat)
}