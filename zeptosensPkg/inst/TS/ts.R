library(zeptosensPkg)
library(zeptosensUtils)
library(paxtoolsr)

nDose=10
nProt=152
maxDist=1 # changing this value requires additional work to compute product(wk). This is not a priority
cellLine="cov318"

#read proteomic response
inputFile <- file.path("inst", "targetScoreData", paste0(cellLine,"_jq1_ave.txt"))
x <- read.table(inputFile, header=TRUE)
#x(dose,prot)
rownames(x) <- x[,1]
x <-x[,-1:-2]
proteomicResponses <- x

targetScoreOutputFile <-"tso.txt"
matrixWkOutputFile <- "wk.txt"
nPerm=30

#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

results <- getTargetScore(nDose, nProt, proteomicResponses, 
                          maxDist=1, nPerm,cellLine, 
                          targetScoreOutputFile="inst/TS/tso.txt", 
                          matrixWkOutputFile="inst/TS/wk.txt",
                          targetScoreQValueFile="inst/TS/q.txt", 
                          targetScoreDoseFile="inst/TS/tsd.txt")
