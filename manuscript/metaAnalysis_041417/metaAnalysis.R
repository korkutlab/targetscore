library(metap)

time <- "24h"
cellline <- "bt474"
dimens <- "3D"
drug <- c("JQ1","iBET151","iBET726","iBET762")
folder <- "../ts_bt474/"

TS1 <- read.table(
  paste(folder,cellline,"_",drug[1],"_",dimens,"_",time,"_TS.txt",sep=""),
  header=F,row.names=1)
P1 <- read.table(
  paste(folder,cellline,"_",drug[1],"_",dimens,"_",time,"_p.txt",sep=""),
  header=T,row.names=1)

TS2 <- read.table(
  paste(folder,cellline,"_",drug[2],"_",dimens,"_",time,"_TS.txt",sep=""),
  header=F,row.names=1)
P2 <- read.table(
  paste(folder,cellline,"_",drug[2],"_",dimens,"_",time,"_p.txt",sep=""),
  header=T,row.names=1)

TS3 <- read.table(
  paste(folder,cellline,"_",drug[3],"_",dimens,"_",time,"_TS.txt",sep=""),
  header=F,row.names=1)
P3 <- read.table(
  paste(folder,cellline,"_",drug[3],"_",dimens,"_",time,"_p.txt",sep=""),
  header=T,row.names=1)

TS4 <- read.table(
  paste(folder,cellline,"_",drug[4],"_",dimens,"_",time,"_TS.txt",sep=""),
  header=F,row.names=1)
P4 <- read.table(
  paste(folder,cellline,"_",drug[4],"_",dimens,"_",time,"_p.txt",sep=""),
  header=T,row.names=1)

colnames(P1) <- "Pval1"
colnames(P2) <- "Pval2"
colnames(P3) <- "Pval3"
colnames(P4) <- "Pval4"

PVex1 <- P1[!is.na(P1$Pval1),]
PVex2 <- P2[!is.na(P2$Pval2),]
PVex3 <- P3[!is.na(P3$Pval3),]
PVex4 <- P4[!is.na(P4$Pval4),]

Pv <- cbind(P1,P2,P3,P4)
#rownames(Pv) <- rownames(P1)
PVex <- cbind(PVex1,PVex2,PVex3,PVex4)
rownames(PVex)<- rownames(P1)[which(!is.na(P1$Pval1))]

CumP <- array(0:0,dim=c(nrow(PVex)))
for (i in 1:nrow(PVex)){
  print(PVex[i,])
  pp <- as.vector(c(PVex[i,1],PVex[i,2],PVex[i,3],PVex[i,4]))
  wght <- as.vector(c(0.5,1,1,1))
  CumP[i] <- sumz(pp,weights=wght)$p 
}
rownames(CumP) <- rownames(PVex)

Q1<- p.adjust(PVex1,method="BH",n=nrow(PVex))
Q2<- p.adjust(PVex2,method="BH",n=nrow(PVex))
Q3<- p.adjust(PVex3,method="BH",n=nrow(PVex))
Q4<- p.adjust(PVex4,method="BH",n=nrow(PVex))
Q<- p.adjust(CumP,method="BH",n=nrow(PVex))

Qall <- cbind(Q1,Q2,Q3,Q4,Q)
TSavg <- (TS1[!is.na(P1$Pval1),1]+
            2*TS2[!is.na(P1$Pval1),1]+
            2*TS3[!is.na(P1$Pval1),1]+
            2*TS4[!is.na(P1$Pval1),1])/7

rownames(Qall)<- rownames(PVex)
QTS <- cbind(TSavg,Q)
rownames(QTS)<- rownames(PVex)
Qall[order(Qall[,5]),]
QTSsorted <- QTS[order(QTS[,1]),]

QTSsortedsig <- QTSsorted[which(QTSsorted[,2] < 0.05),]
QTSsortedsig2 <- cbind(rownames(QTSsortedsig),QTSsortedsig[,1],QTSsortedsig[,1],QTSsortedsig[,2])
#Qall[sort]
colnames(QTSsortedsig2) <- c("ID1","changeInclude","TSavg","Meta_FDRp")

write.table(QTSsortedsig2, file=paste("Meta_",cellline,"_",time,"_",dimens,".txt",sep=""),
            quote=F,row.names=F)
rm(list=ls())
