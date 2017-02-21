#generate sif files for selected 
wk <- read.table("wks.txt",header=T,row.names = 1)
prot_names <- read.table("prot_names.txt",header=T,row.names = 1)
rownames(wk) <- row.names(prot_names)
colnames(wk) <- row.names(prot_names)
inter <-which(wk!=0,arr.ind=T) 
sifi <- array(0,dim=c(nrow(inter),3))
for (i in 1:nrow(inter)){
    sifi[i,1]<- rownames(wk[inter[i,][1],][inter[i,][2]])
    sifi[i,2] <- as.numeric(wk[inter[i,][1],][inter[i,][2]])
    sifi[i,3] <- colnames(wk[inter[i,][1],][inter[i,][2]])
    }

write.table(file="sif_all.txt",sifi,quote=F,col.names=F,row.names=F)
rm(list = ls())
