#generate sif files for selected 
wk <- read.table("wks.txt",header=T,row.names = 1)
prot_names <- read.table("protlist1.txt",header=T,row.names = 1)
rownames(wk) <- row.names(prot_names)
colnames(wk) <- row.names(prot_names)
upst <- matrix(wk[,colnames(wk) %in% "STAT3pY705"])
#upst <- as.matrix(wk)
sifu<- upst[subset(upst!=0),]
#inter <-subset(upst!=0) 
sifd <- rep("STAT3pY705",each=nrow(sifu))
sif <- cbind(rownames(wk[which(upst!=0),]),upst[subset(upst!=0),],"STAT3pY705")
write.table(file="sif_STAT3pY705.txt",sif,quote=F,col.names=F,row.names = F)
rm(list = ls())
