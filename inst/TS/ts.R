library(zeptosensPkg)
library(paxtoolsr)
nDose=10
nProt=152
maxDist=1 # changing this value requires additional work to compute product(wk). This is not a priority
cellLine="cov318"

#read proteomic response
inputFile <- file.path("~/zeptosenspkg/inst", "targetScoreData", paste0(cellLine,"_jq1_ave.txt"))
x <- read.table(inputFile, header=TRUE)
#x(dose,prot)
rownames(x) <- x[,1]
x <-x[,-1:-2]
proteomicResponses <- x

targetScoreOutputFile <-"tso.txt"
matrixWkOutputFile <- "wk.txt"
nPerm=2
results <- getTargetScore(nDose, nProt, proteomicResponses, 
                          maxDist=1, nPerm,cellLine, 
                          targetScoreOutputFile="~/zeptosenspkg/inst/TS/tso.txt", 
                          matrixWkOutputFile="~/zeptosenspkg/inst/TS/wk.txt",
                          targetScoreQValueFile="~/zeptosenspkg/inst/TS/q.txt", 
                          targetScoreDoseFile="~/zeptosenspkg/inst/TS/tsd.txt")
results