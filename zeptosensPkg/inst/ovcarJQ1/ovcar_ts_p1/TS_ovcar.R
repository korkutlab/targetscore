library(xlsx)
library(zeptosensUtils)
library(zeptosensPkg)


######
nDose=10
nProt=152
maxDist=1 # changing this value requires additional work to compute product(wk). This is not a priority
nSample=4
cellLine1="ovcar3"
cellLine2="ovcar4"
cellLine3="igrov1"
cellLine4="cov318"

#read proteomic response for cellline1
inputFile <- file.path("../../../inst", "ovcarJQ1/data", paste0(cellLine1,"_jq1_ave.txt"))
x_1 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_1) <- x_1[,1]
x_1 <- x_1[,-1:-2]

#read proteomic response for cellline2
inputFile <- file.path("../../../inst", "ovcarJQ1/data", paste0(cellLine2,"_jq1_ave.txt"))
x_2 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_2) <- x_2[,1]
x_2 <- x_2[,-1:-2]

#read proteomic response for cellline3
inputFile <- file.path("../../../inst", "ovcarJQ1/data", paste0(cellLine3,"_jq1_ave.txt"))
x_3 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_3) <- x_3[,1]
x_3 <- x_3[,-1:-2]

#read proteomic response for cellline4
inputFile <- file.path("../../../inst", "ovcarJQ1/data", paste0(cellLine4,"_jq1_ave.txt"))
x_4 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_4) <- x_4[,1]
x_4 <- x_4[,-1:-2]


nX <- rbind(x_1,x_2,x_3,x_4)


stdev <- sampSdev(nSample=nSample,nProt=nProt,nDose=nDose,nX=nX)



#TS for cellLine1

targetScoreOutputFile <-paste0(cellLine1,"_TS.txt")
matrixWkOutputFile <- "wk_1.txt"
signedMatrixWkOutputFile <- "wks.txt"
nPerm=3
maxDist <- 1
proteomicResponses_1 <- x_1
for(i in 1:nProt){
    for (j in 1:nDose){
        proteomicResponses_1[j,i] <- (x_1[j,i]/stdev[i])      
    }
}

length(proteomicResponses_1)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses_1, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine1, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine1,"_q.txt"), 
                          targetScoreDoseFile=paste0(cellLine1,"_TS_d.txt"),
                          targetScorePValueFile=paste0(cellLine1,"_p.txt"),
                          verbose=FALSE,fsFile="fs.txt",
                          signedMatrixWkOutputFile=signedMatrixWkOutputFile)


#TS for cellLine2

targetScoreOutputFile <-paste0(cellLine2,"_TS.txt")
matrixWkOutputFile <- "wk_2.txt"
nPerm=30
maxDist <- 1
proteomicResponses_2 <- x_2
for(i in 1:nProt){
    for (j in 1:nDose){
        proteomicResponses_2[j,i] <- (x_2[j,i]/stdev[i])      
    }
}

length(proteomicResponses_2)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses_2, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine2, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine2,"_q.txt"), 
                          targetScoreDoseFile=paste0(cellLine2,"_TS_d.txt"),
                          targetScorePValueFile=paste0(cellLine2,"_p.txt"),
                          verbose=FALSE,fsFile="fs.txt")



#TS for cellLine3

targetScoreOutputFile <-paste0(cellLine3,"_TS.txt")
matrixWkOutputFile <- "wk_3.txt"
nPerm=30
maxDist <- 1
proteomicResponses_3 <- x_3
for(i in 1:nProt){
    for (j in 1:nDose){
        proteomicResponses_3[j,i] <- (x_3[j,i]/stdev[i])      
    }
}

length(proteomicResponses_3)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses_3, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine3, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine3,"_q.txt"), 
                          targetScoreDoseFile=paste0(cellLine3,"_TS_d.txt"),
                          targetScorePValueFile=paste0(cellLine3,"_p.txt"),
                          verbose=FALSE,fsFile="fs.txt")



#TS for CellLine4


targetScoreOutputFile <-paste0(cellLine4,"_TS.txt")
matrixWkOutputFile <- "wk_4.txt"
nPerm=30
maxDist <- 1
proteomicResponses_4 <- x_4
for(i in 1:nProt){
    for (j in 1:nDose){
        proteomicResponses_4[j,i] <- (x_4[j,i]/stdev[i])      
    }
}

length(proteomicResponses_4)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses_4, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine4, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine4,"_q.txt"), 
                          targetScoreDoseFile=paste0(cellLine4,"_TS_d.txt"),
                          targetScorePValueFile=paste0(cellLine4,"_p.txt"),
                          verbose=FALSE,fsFile="fs.txt")


rm(list = ls())
