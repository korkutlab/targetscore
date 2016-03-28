#' Read Zeptosens Export File 
#' 
#' @param data Zeptosens data read in by readZeptosensXls
#' @param sampleNameEntries FIXME, chracter string, perturbation conditions, #cell line etc (sampleNameEntries <- c('sNumber','sample', 'treatment', NA,'time',NA,NA,NA))#
#' 
#' @details 
#' 
#' Format: sNumber_sample_treatment_dose_time_replicate_date_notes
#' 
#' @examples 
#' tmp <- readZeptosensExport(tmp, c('sNumber', NA, NA, 'sample', 'treatment', 'time', 'replicate'))
#' 
#' @concept zeptosens
#' @export
readZeptosensExport <- function(data, sampleNameEntries = NULL) {
    message("NOTE: Sample Name in *_RFI_Export_Table.xls is assumed to be in the format sNumber_sample_treatment_dose_time_replicate_date_notes. Use 'samplenNameEntries' parameter to modify as described in documentation.")
    
    splitSampleNames <- splitSampleNames(data[, "Sample.Name"], sampleNameEntries)
    
    df <- NULL
    
    allCols <- c("sNumber", "sample", "treatment", "dose", "time", "replicate", "date", 
        "notes", "antibody", "readout", "cv", "quality")
    
    for (i in 1:nrow(data)) {
        tmpRow <- NULL
        
        for (col in allCols) {
            if (col == "sNumber") {
                tmpRow <- c(tmpRow, sNumber = splitSampleNames[i, col])
            }
            
            if (col == "sample") {
                tmpRow <- c(tmpRow, sample = splitSampleNames[i, col])
            }
            
            if (col == "treatment") {
                
                tmpRow <- c(tmpRow, treatment = splitSampleNames[i, col])
            }
            
            if (col == "dose") {
                tmpRow <- c(tmpRow, dose = splitSampleNames[i, col])
            }
            
            if (col == "time") {
                
                tmpRow <- c(tmpRow, time = splitSampleNames[i, col])
            }
            
            if (col == "replicate") {
                tmpRow <- c(tmpRow, replicate = splitSampleNames[i, col])
            }
            
            if (col == "notes") {
                tmpRow <- c(tmpRow, notes = splitSampleNames[i, col])
            }
            
            if (col == "antibody") {
                tmpRow <- c(tmpRow, antibody = data[i, "Analyte"])
            }
            
            if (col == "readout") {
                tmpRow <- c(tmpRow, readout = data[i, "RFI"])
            }
            
            if (col == "cv") {
                tmpRow <- c(tmpRow, cv = data[i, "RFI.CV"])
            }
            
            if (col == "quality") {
                tmpRow <- c(tmpRow, quality = data[i, "Class."])
            }
        }
        
        df <- rbind(df, tmpRow)
    }
    
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    # Convert appropriate columns to numeric
    numCols <- c("readout", "cv")
    df[numCols] <- sapply(df[numCols], as.numeric)
    
    rownames(df) <- 1:nrow(df)
    
    return(df)
}

# #' Read Zeptosens Export File #' #' @param inputFile TBA #' @param
# sampleNameEntries TBA #' #' @examples #' tmp <-
# readZeptosensExport('inst/dataInst/R006_RFI_Export_Table.xls',
# c('sampleNameNumber', 'NA', 'NA', 'sampleName', 'treatment', 'time',
# 'replicate')) #' #' @concept zeptosensPkg readZeptosensExportDEPRECATED <-
# function(data, sampleNameEntries=NULL) { message('NOTE: Sample Name in
# *_RFI_Export_Table.xls is assumed to be in the format
# sNumber_sampleName_drug_time_dose_replicateNumber. Use 'samplenNameEntries'
# parameter to modify as described in documentation.') splitSampleNames <-
# strsplit(data[, 'Sample.Name'], '_') df <- NULL for(i in 1:nrow(data)) {
# sampleNameNumberIdx <- which(sampleNameEntries == 'sampleNameNumber')
# sampleNameIdx <- which(sampleNameEntries == 'sampleName') treatmentIdx <-
# which(sampleNameEntries == 'treatment') timeIdx <- which(sampleNameEntries ==
# 'time') replicateIdx <- which(sampleNameEntries == 'replicate') # TODO: More
# generalized tmpDf <- data.frame( sNumber=data[i, 'Sample.'],
# sampleNameNumber=splitSampleNames[[i]][sampleNameNumberIdx],
# sampleName=splitSampleNames[[i]][sampleNameIdx],
# treatment=splitSampleNames[[i]][treatmentIdx],
# time=splitSampleNames[[i]][timeIdx],
# replicate=splitSampleNames[[i]][replicateIdx], antibody=data[i, 'Analyte'],
# readout=data[i, 'RFI'], cv=data[i, 'RFI.CV'], quality=data[i, 'Class.'],
# stringsAsFactors=FALSE ) df <- rbind(df, tmpDf) } return(df) } #' Read
# Intermediate Export File #' #' @param nrep number of replicates #' @param ncond
# number of conditions #' @param ntime number of timepoints #' @param nab number
# of antibodies #' #' @expamples #' #TODO #' #' @concept zeptosensPkg #
# readIntermediateExport <- function(data, outputFile=NULL, nrep=NULL,
# ncond=NULL, ntime=NULL, nab=NULL) { options('digits'=3) nrep=2 ncond=3 ntime=4
# nab=72 ctr=ncond*ntime*nrep ct=ncond*ntime #data <- read.table('a2058.txt')
# annot <- read.table('ab_index.txt') antibodyName geneName
# postTranslationModification activityState normal <- read.table('cellmix.txt')
# for (i in 1:length(data[,1])) { if(data[i, 5] == 'Poor') { data[i, 3] = 0
# data[i, 4] = 0 } } data.matrix <- as.matrix(data) annot.matrix <-
# as.matrix(annot) datar <- array(dim=c(nab*ncond*ntime,7)) for(ab in 1:nab){ for
# (k in 1:ct) { # print(k) datar[ct*(ab-1)+k,1]=data.matrix[ctr*(ab-1)+k,2]
# datar[ct*(ab-1)+k,2]=annot.matrix[ab,1] datar[ct*(ab-1)+k,3]=annot.matrix[ab,2]
# datar[ct*(ab-1)+k,4]=annot.matrix[ab,3] datar[ct*(ab-1)+k,5]=annot.matrix[ab,4]
# # datar[ct*(ab-1)+k,5]=data.matrix[ctr*(ab-1)+k,2] #
# datar[ct*(ab-1)+k,6]=data.matrix[ctr*(ab-1)+k,3] #
# datar[ct*(ab-1)+k,7]=data.matrix[ctr*(ab-1)+ct+k,3]
# datar[ct*(ab-1)+k,6]=data[ctr*(ab-1)+k,4]#/(normal[ab,4])
# datar[ct*(ab-1)+k,7]=data[ctr*(ab-1)+ct+k,4]#/(normal[ab,4]) } }
# write.table(datar, file=outputFile, col.names=TRUE, row.names=FALSE, sep='\t',
# quote=FALSE) } #' Plot Zeptosens Replicates #' #' @param outputFile #' @param
# plotZeptosensReplicates <- function(data, outputFile='replicates.pdf', ) {
# g_datar=datar[which(as.numeric(datar[,6]) > 0 & as.numeric(datar[,7]) > 0 ),]
# cc <- cor(as.numeric(g_datar[,6]),as.numeric(g_datar[,7]))
# pdf('replicates.pdf') plot(g_datar[,6],g_datar[,7],xlab='replicate
# 1',ylab='replicate 2') title(main=paste('a2058 replicates \nCC=',cc))
# dev.off() pdf('replicates_2.pdf') plot(g_datar[,6],g_datar[,7],xlab='replicate
# 1',ylab='replicate 2',xlim=c(0,2),ylim=c(0,2)) title(main=paste('a2058
# replicates\nCC=',cc)) dev.off() rm(list=ls()) } 
