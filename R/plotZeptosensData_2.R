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
plotZeptosensData <- function(normArray1, normArray2, antibodies, excludeTreatments, 
                              ctrlTreatment, ctrlTime, plotColNames, plotRowNames, plotColors, plotDir) {
    allDrugs <- unique(normArray1[, "treatment"])
    tmpDrugs <- allDrugs[!(unique(normArray1[, "treatment"]) %in% excludeTreatments)]
    
    curDrug <- tmpDrugs[1]
    
    for(curAntibody in antibodies) {
#        plotRows <- ceiling(sqrt(length(antibodies)))
#        plotCols <- ceiling(sqrt(length(antibodies)))
        
        pdf(paste0(plotDir, curAntibody, ".pdf"), width=20, height=8.5)
        
#        par(mfrow=c(plotRows, plotCols))
        
        par(mar=c(2,2,1.5,1), oma=c(0,0,2,0))
#        for(curAntibody in antibodies) {
            tmpArray1Entries <- normArray1[which(normArray1[, "antibody"] == curAntibody), ]
#            tmpArray2Entries <- normArray2[which(normArray2[, "antibody"] == curAntibody), ]
            
            idx <- intersect(which(tmpArray1Entries[, "antibody"] == curAntibody), which(tmpArray1Entries[, "treatment"] == curDrug))
            
            t0Idx <- intersect(which(tmpArray1Entries[, "treatment"] == ctrlTreatment), which(tmpArray1Entries[, "time"] == ctrlTime))
            allIdx <- c(t0Idx, idx)
            
#            entriesWithCtrl <- cbind(tmpArray1Entries[allIdx, "normReadout"], tmpArray2Entries[allIdx, "normReadout"])
             entriesWithCtrl <- tmpArray1Entries[allIdx, "normReadout"]
            for(curDrug in tmpDrugs) {
                colnames(entriesWithCtrl) <- plotColNames
                rownames(entriesWithCtrl) <- plotRowNames
                tmpMat <- t(entriesWithCtrl)
                
                plot(tmpMat, main=curAntibody, beside=TRUE, col=plotColors, cex.main=0.9, cex.axis=0.7)                
            }
 
#        }
        title(main=curDrug, outer=TRUE)
        
        dev.off()
    }
    
}