#generate sif files for selected 
wk <- read.table("wks.txt",header=T,row.names = 1)
prot_names <- read.table("prot_names.txt",header=T,row.names = 1)
rownames(wk) <- row.names(prot_names)
colnames(wk) <- row.names(prot_names)
upst <- matrix(wk[,colnames(wk) %in% "PKCalpha"])
colnames(upst)<- "PKCalpha"
rownames(upst) <- as.character(row.names(prot_names))
sifu<- upst[subset(upst!=0),,drop=FALSE]
sifd <- rep("PKCalpha",each=nrow(sifu))
sif <- cbind(sifu,sifd)
write.table(file="sif_PKCalpha.txt",sif,quote=F,col.names=F)
rm(list = ls())
