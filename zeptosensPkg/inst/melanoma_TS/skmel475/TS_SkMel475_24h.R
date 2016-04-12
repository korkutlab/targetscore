library(xlsx)
library(zeptosensUtils)
library(zeptosensPkg)

## Normalization

inputFile <- file.path("inst", "melanoma_TS/skmel475", "R009_RFI_Export_Table_skmel118_skmel475.xlsx")
#inputFile <- "R009_RFI_Export_Table_skmel118_skmel475.xlsx"
tmpDat <- readZeptosensXls(inputFile)


# Example Sample Name: S001_cell line_Melanoma_A2058 _DMSO_1hr_rep1

sampleNameEntries <- c("sNumber", "sample", "treatment", "time", "replicate", "date", NA)
dat <- readZeptosensExport(tmpDat, sampleNameEntries)


colClasses <- c("character","character","numeric")
inputFile <- file.path("inst", "melanoma_TS/skmel475", "R009_total_protein_level_skmel118_skmel475.xlsx")

tplDat <- read.xlsx2(inputFile , colClasses=colClasses, stringsAsFactors=F, sheetIndex=1)

# This is from the TPL file
nSamples <- 64

antibodyNum <- length(unique(dat[, "antibody"]))

mel_ttl_prot_norm <- totalProteinNormalization(dat, tplDat, nSamples, antibodyNum)

mel475 <- mel_ttl_prot_norm[which(grepl("Mel-475", mel_ttl_prot_norm[, "sample"], ignore.case=TRUE)),]

mel118 <- mel_ttl_prot_norm[which(grepl("Mel-118", mel_ttl_prot_norm[, "sample"], ignore.case=TRUE)),]

mel475_24 <- mel475[which(mel475$time == "24hr"),]
nSamples <- nrow(mel475_24[which(mel475_24$antibody == "Akt"),])

array1 <- mel475_24
normalFactor <- "DMSO"

mel475_ttl_unprtbed_norm <- unperturbedNormalization(array1, normalFactor, nSamples, antibodyNum)

# log2 transformation: 
#mel475_ttl_unprtbed_norm$readout <- log2(mel475_ttl_unprtbed_norm$readout)

# Treatment 901 is the MEKi 

#mel475_ttl_unprtbed_norm

repCnt <- 3 

which(mel475_ttl_unprtbed_norm$replicate == "rep1" & mel475_ttl_unprtbed_norm$treatment != "DMSO")

rep1 <- mel475_ttl_unprtbed_norm[mel475_ttl_unprtbed_norm$replicate == "rep1" & mel475_ttl_unprtbed_norm$treatment != "DMSO", "readout"]
rep2 <- mel475_ttl_unprtbed_norm[mel475_ttl_unprtbed_norm$replicate == "rep2" & mel475_ttl_unprtbed_norm$treatment != "DMSO", "readout"]
rep3 <- mel475_ttl_unprtbed_norm[mel475_ttl_unprtbed_norm$replicate == "rep3" & mel475_ttl_unprtbed_norm$treatment != "DMSO", "readout"]

y1 <- rbind(rep1, rep2, rep3)

y2 <- mel475_ttl_unprtbed_norm[mel475_ttl_unprtbed_norm$replicate == "rep3" & mel475_ttl_unprtbed_norm$treatment != "DMSO", ]
avg <- colSums(y1) / repCnt

y2$readout <- log2(avg)
y2$replicate <- rep("avg", nrow(y2))


y3 <- matrix(y2$readout, 1, antibodyNum)

rownames(y3) <- "MEKi24_Skmel475"
colnames(y3) <- unique(y2$antibody)
inputFile <- file.path("inst", "melanoma_TS/skmel475", "skmel475_meki_ave_24.txt")
write.table(y3,file=inputFile,quote=F,sep= "\t")
test <- read.delim("inst/targetscoreData/antibodyMap.txt")


######
nDose=1
nProt=71
maxDist=1 # changing this value requires additional work to compute product(wk). This is not a priority
cellLine="skmel475"

#read proteomic response
inputFile <- file.path("inst", "melanoma_TS/skmel475", paste0(cellLine,"_meki_ave_24.txt"))
x <- read.table(inputFile, header=TRUE,sep="\t",check.names=FALSE)
#x(dose,prot)
#rownames(x) <- x[,1]
#x <-x[,-1]
x <- x[,-72]

proteomicResponses <- x

targetScoreOutputFile <-"inst/melanoma_TS/skmel475/tso_24.txt"
matrixWkOutputFile <- "inst/melanoma_TS/skmel475/wk.txt"
nPerm=3
maxDist <- 1
length(proteomicResponses)
#results <- calcTargetScore(nDose, nProt, proteomicResponses, maxDist = 1, cellLine)

set.seed(1)

results <- getTargetScore(nDose=nDose, 
                          nProt=nProt, 
                          proteomicResponses=proteomicResponses, 
                          maxDist=maxDist, 
                          nPerm=nPerm,
                          cellLine=cellLine, 
                          targetScoreOutputFile=targetScoreOutputFile, 
                          matrixWkOutputFile=matrixWkOutputFile,
                          targetScoreQValueFile="inst/melanoma_TS/skmel475/q_24.txt", 
                          targetScoreDoseFile="inst/melanoma_TS/skmel475/tsd2_24.txt",
                          verbose=TRUE,
                          tsFactor=1,
                          fsFile="inst/targetScoreData/fs_mskcc.txt")
