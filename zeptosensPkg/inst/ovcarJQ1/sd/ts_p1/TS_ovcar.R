library(xlsx)
library(zeptosensUtils)
library(zeptosensPkg)
source('~/zeptosenspkg/zeptosensPkg/R/sampSdev.R')

######
nDose=10
nProt=152
maxDist=1 # changing this value requires additional work to compute product(wk). This is not a priority
nSample=12
cellLine1="ovcar3"
cellLine2="ovcar4"
cellLine3="igrov1"
cellLine4="cov318"

#read proteomic response for cellline1
inputFile <- file.path("../../../../inst", "ovcarJQ1/sd/data", paste0(cellLine1,"_jq1_ci.txt"))
x_1 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_1) <- x_1[,1]
x_1 <- x_1[,-1:-2]
x1_1 <- x_1[1:10,]
x1_2 <- x_1[11:20,]
x1_3 <- x_1[21:30,]
#read proteomic response for cellline2
inputFile <- file.path("../../../../inst", "ovcarJQ1/sd/data", paste0(cellLine2,"_jq1_ci.txt"))
x_2 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_2) <- x_2[,1]
x_2 <- x_2[,-1:-2]
x2_1 <- x_2[1:10,]
x2_2 <- x_2[11:20,]
x2_3 <- x_2[21:30,]
#read proteomic response for cellline3
inputFile <- file.path("../../../../inst", "ovcarJQ1/sd/data", paste0(cellLine3,"_jq1_ci.txt"))
x_3 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_3) <- x_3[,1]
x_3 <- x_3[,-1:-2]
x3_1 <- x_2[1:10,]
x3_2 <- x_2[11:20,]
x3_3 <- x_2[21:30,]
#read proteomic response for cellline4
inputFile <- file.path("../../../../inst", "ovcarJQ1/sd/data", paste0(cellLine4,"_jq1_ci.txt"))
x_4 <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
rownames(x_4) <- x_4[,1]
x_4 <- x_4[,-1:-2]
x4_1 <- x_4[1:10,]
x4_2 <- x_4[11:20,]
x4_3 <- x_4[21:30,]

nX <- rbind(x1_1,x1_2,x1_3,x2_1,x2_2,x2_3,x3_1,x3_2,x3_3,x4_1,x4_2,x4_3)


stdev <- sampSdev(nSample=nSample,nProt=nProt,nDose=nDose,nX=nX)



#TS for cellLine1
#rep_1
targetScoreOutputFile <-paste0(cellLine1,"_TS_1.txt")
matrixWkOutputFile <- "wk_1_1.txt"
nPerm=0
maxDist <- 1
proteomicResponses1_1 <- x1_1
for(i in 1:nProt){
    for (j in 1:nDose){
        proteomicResponses1_1[j,i] <- (x1_1[j,i]/stdev[i])      
    }
}

length(proteomicResponses1_1)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses1_1, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine1, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine1,"_q_1.txt"), 
                          targetScoreDoseFile=paste0(cellLine1,"_TS_d_1.txt"),
                          targetScorePValueFile=paste0(cellLine1,"_p_1.txt"),
                          verbose=FALSE,fsFile="fs.txt")

#####
#rep2
targetScoreOutputFile <-paste0(cellLine1,"_TS_2.txt")
matrixWkOutputFile <- "wk1_2.txt"
nPerm=0
maxDist <- 1
proteomicResponses1_2 <- x1_2
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses1_2[j,i] <- (x1_2[j,i]/stdev[i])      
  }
}

length(proteomicResponses1_2)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses1_2, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine1, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine1,"_q_2.txt"), 
                          targetScoreDoseFile=paste0(cellLine1,"_TS_d_2.txt"),
                          targetScorePValueFile=paste0(cellLine1,"_p_2.txt"),
                          verbose=FALSE,fsFile="fs.txt")

###
#rep_3
targetScoreOutputFile <-paste0(cellLine1,"_TS_3.txt")
matrixWkOutputFile <- "wk1_3.txt"
nPerm=0
maxDist <- 1
proteomicResponses1_3 <- x1_3
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses1_3[j,i] <- (x1_3[j,i]/stdev[i])      
  }
}

length(proteomicResponses1_3)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses1_3, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine1, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine1,"_q_3.txt"), 
                          targetScoreDoseFile=paste0(cellLine1,"_TS_d_3.txt"),
                          targetScorePValueFile=paste0(cellLine1,"_p_3.txt"),
                          verbose=FALSE,fsFile="fs.txt")
tsd1_1 <- read.table(paste(cellLine1,"_TS_d_1.txt",sep=""),header=T)
tsd1_2 <- read.table(paste(cellLine1,"_TS_d_2.txt",sep=""),header=T)
tsd1_3 <- read.table(paste(cellLine1,"_TS_d_3.txt",sep=""),header=T)

ts_sd_c1 <- 100*tsd1_1/tsd1_1

for(i in 1:nDose){
  for (j in 1:nProt){
    ts_sd_c1[i,j]=sd(c(tsd1_1[i,j],tsd1_2[i,j],tsd1_3[i,j]))
  }
}
#TS for cellLine2
#rep_1
targetScoreOutputFile <-paste0(cellLine2,"_TS_1.txt")
matrixWkOutputFile <- "wk_2_1.txt"
nPerm=0
maxDist <- 1
proteomicResponses2_1 <- x2_1
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses2_1[j,i] <- (x2_1[j,i]/stdev[i])      
  }
}

length(proteomicResponses2_1)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses2_1, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine2, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine2,"_q_1.txt"), 
                          targetScoreDoseFile=paste0(cellLine2,"_TS_d_1.txt"),
                          targetScorePValueFile=paste0(cellLine2,"_p_1.txt"),
                          verbose=FALSE,fsFile="fs.txt")

#####
#rep2
targetScoreOutputFile <-paste0(cellLine2,"_TS_2.txt")
matrixWkOutputFile <- "wk2_2.txt"
nPerm=0
maxDist <- 1
proteomicResponses2_2 <- x2_2
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses2_2[j,i] <- (x2_2[j,i]/stdev[i])      
  }
}

length(proteomicResponses2_2)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses2_2, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine2, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine2,"_q_2.txt"), 
                          targetScoreDoseFile=paste0(cellLine2,"_TS_d_2.txt"),
                          targetScorePValueFile=paste0(cellLine2,"_p_2.txt"),
                          verbose=FALSE,fsFile="fs.txt")

###
#rep_3
targetScoreOutputFile <-paste0(cellLine2,"_TS_3.txt")
matrixWkOutputFile <- "wk2_3.txt"
nPerm=0
maxDist <- 1
proteomicResponses2_3 <- x2_3
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses2_3[j,i] <- (x2_3[j,i]/stdev[i])      
  }
}

length(proteomicResponses2_3)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses2_3, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine2, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellLine2,"_q_3.txt"), 
                          targetScoreDoseFile=paste0(cellLine2,"_TS_d_3.txt"),
                          targetScorePValueFile=paste0(cellLine2,"_p_3.txt"),
                          verbose=FALSE,fsFile="fs.txt")
tsd2_1 <- read.table(paste(cellLine2,"_TS_d_1.txt",sep=""),header=T)
tsd2_2 <- read.table(paste(cellLine2,"_TS_d_2.txt",sep=""),header=T)
tsd2_3 <- read.table(paste(cellLine2,"_TS_d_3.txt",sep=""),header=T)

ts_sd_c2 <- 100*tsd2_1/tsd2_1

for(i in 1:nDose){
  for (j in 1:nProt){
    ts_sd_c2[i,j]=sd(c(tsd2_1[i,j],tsd2_2[i,j],tsd2_3[i,j]))
  }
}


#TS for cellline3
#rep_1
targetScoreOutputFile <-paste0(cellline3,"_TS_1.txt")
matrixWkOutputFile <- "wk_3_1.txt"
nPerm=0
maxDist <- 1
proteomicResponses3_1 <- x3_1
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses3_1[j,i] <- (x3_1[j,i]/stdev[i])      
  }
}

length(proteomicResponses3_1)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses3_1, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellline3, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellline3,"_q_1.txt"), 
                          targetScoreDoseFile=paste0(cellline3,"_TS_d_1.txt"),
                          targetScorePValueFile=paste0(cellline3,"_p_1.txt"),
                          verbose=FALSE,fsFile="fs.txt")

#####
#rep2
targetScoreOutputFile <-paste0(cellline3,"_TS_2.txt")
matrixWkOutputFile <- "wk3_2.txt"
nPerm=0
maxDist <- 1
proteomicResponses3_2 <- x3_2
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses3_2[j,i] <- (x3_2[j,i]/stdev[i])      
  }
}

length(proteomicResponses3_2)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses3_2, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellline3, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellline3,"_q_2.txt"), 
                          targetScoreDoseFile=paste0(cellline3,"_TS_d_2.txt"),
                          targetScorePValueFile=paste0(cellline3,"_p_2.txt"),
                          verbose=FALSE,fsFile="fs.txt")

###
#rep_3
targetScoreOutputFile <-paste0(cellline3,"_TS_3.txt")
matrixWkOutputFile <- "wk3_3.txt"
nPerm=0
maxDist <- 1
proteomicResponses3_3 <- x3_3
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses3_3[j,i] <- (x3_3[j,i]/stdev[i])      
  }
}

length(proteomicResponses3_3)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses3_3, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellline3, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellline3,"_q_3.txt"), 
                          targetScoreDoseFile=paste0(cellline3,"_TS_d_3.txt"),
                          targetScorePValueFile=paste0(cellline3,"_p_3.txt"),
                          verbose=FALSE,fsFile="fs.txt")
tsd3_1 <- read.table(paste(cellline3,"_TS_d_1.txt",sep=""),header=T)
tsd3_2 <- read.table(paste(cellline3,"_TS_d_2.txt",sep=""),header=T)
tsd3_3 <- read.table(paste(cellline3,"_TS_d_3.txt",sep=""),header=T)

ts_sd_c3 <- 100*tsd3_1/tsd3_1

for(i in 1:nDose){
  for (j in 1:nProt){
    ts_sd_c3[i,j]=sd(c(tsd3_1[i,j],tsd3_2[i,j],tsd3_3[i,j]))
  }
}
#TS for cellline4
#rep_1
targetScoreOutputFile <-paste0(cellline4,"_TS_1.txt")
matrixWkOutputFile <- "wk_4_1.txt"
nPerm=0
maxDist <- 1
proteomicResponses4_1 <- x4_1
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses4_1[j,i] <- (x4_1[j,i]/stdev[i])      
  }
}

length(proteomicResponses4_1)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses4_1, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellline4, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellline4,"_q_1.txt"), 
                          targetScoreDoseFile=paste0(cellline4,"_TS_d_1.txt"),
                          targetScorePValueFile=paste0(cellline4,"_p_1.txt"),
                          verbose=FALSE,fsFile="fs.txt")

#####
#rep2
targetScoreOutputFile <-paste0(cellline4,"_TS_2.txt")
matrixWkOutputFile <- "wk4_2.txt"
nPerm=0
maxDist <- 1
proteomicResponses4_2 <- x4_2
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses4_2[j,i] <- (x4_2[j,i]/stdev[i])      
  }
}

length(proteomicResponses4_2)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses4_2, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellline4, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellline4,"_q_2.txt"), 
                          targetScoreDoseFile=paste0(cellline4,"_TS_d_2.txt"),
                          targetScorePValueFile=paste0(cellline4,"_p_2.txt"),
                          verbose=FALSE,fsFile="fs.txt")

###
#rep_3
targetScoreOutputFile <-paste0(cellline4,"_TS_3.txt")
matrixWkOutputFile <- "wk4_3.txt"
nPerm=0
maxDist <- 1
proteomicResponses4_3 <- x4_3
for(i in 1:nProt){
  for (j in 1:nDose){
    proteomicResponses4_3[j,i] <- (x4_3[j,i]/stdev[i])      
  }
}

length(proteomicResponses4_3)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses4_3, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellline4, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile=paste0(cellline4,"_q_3.txt"), 
                          targetScoreDoseFile=paste0(cellline4,"_TS_d_3.txt"),
                          targetScorePValueFile=paste0(cellline4,"_p_3.txt"),
                          verbose=FALSE,fsFile="fs.txt")
tsd4_1 <- read.table(paste(cellline4,"_TS_d_1.txt",sep=""),header=T)
tsd4_2 <- read.table(paste(cellline4,"_TS_d_2.txt",sep=""),header=T)
tsd4_3 <- read.table(paste(cellline4,"_TS_d_3.txt",sep=""),header=T)

ts_sd_c4 <- 100*tsd4_1/tsd4_1

for(i in 1:nDose){
  for (j in 1:nProt){
    ts_sd_c4[i,j]=sd(c(tsd4_1[i,j],tsd4_2[i,j],tsd4_3[i,j]))
  }
}

rm(list = ls())
