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

targetScoreOutputFile <-"inst/melanoma_TS/skmel475/tso_6.txt"
matrixWkOutputFile <- "inst/melanoma_TS/skmel475/wk.txt"
nPerm=30
maxDist <- 1
length(proteomicResponses)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile="inst/melanoma_TS/skmel475/q_6.txt", 
                          targetScoreDoseFile="inst/melanoma_TS/skmel475/tsd2_6.txt",
                          verbose=TRUE,TSfactor=1)
