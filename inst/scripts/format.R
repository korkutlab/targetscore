#script to generate the formatted zeptosens format
nrep=3
ncond=3
ntime=4
nsamp=2
samplename=c("skmel475","skmel118")
normfile="cellmix.txt"

nab=72

filename <-paste(samplename[1],".txt",sep="")
data <- read.table(filename,header=FALSE)
#print(nrow(data))
annot <- read.table("ab_table.txt")
normal <- read.table(normfile)
totalprot <- read.table("total_prot_skmel475.txt", header=T)
snctr<-data.frame(matrix(unlist(strsplit(as.character(data[,8]),"_")),nrow=nrow(data),byrow=T))


#print(annot)
annot_m <- array(dim=c(nrow(data),4))
head <-(c("sample_name","condition","time","replicate","Ab_name","gene_name","posttrans","activity","readout","CV","Quality"))
for (i in 1:nrow(data)){
#	print((as.character(data[i,9])))
#	print((annot))
#    print(((annot[(which(annot[,1]==as.character(data[i,9]))),])))
	annot_m[i,]=as.matrix((annot[(which(annot[,1]==as.character(data[i,9]))),]))
}
#print(noquote(annot_m))
data_formatted <- cbind(snctr[,-1],annot_m,data[,14],data[,15],data[,24])
print(data_formatted)
#print(annot[,1])
#	print((which(data[,9]==as.character(annot[,1]))))
#print(grep(as.character(data[1,9]),annot[,1]))
#print(as.character(data[,9]))
data_r1 <- data_formatted[grep("A",data_formatted[,4]),]
data_r2 <- data_formatted[grep("B",data_formatted[,4]),]
data_r3 <- data_formatted[grep("C",data_formatted[,4]),]
write.table(t(head),file="skmel475_r1.txt",col.names=F,row.names=F,quote=F)
write.table(data_r1,file="skmel475_r1.txt",col.names=F,row.names=F,quote=F,append=T)
write.table(t(head),file="skmel475_r2.txt",col.names=F,row.names=F,quote=F)
write.table(data_r2,file="skmel475_r2.txt",col.names=F,row.names=F,quote=F,append=T)
write.table(t(head),file="skmel475_r3.txt",col.names=F,row.names=F,quote=F)
write.table(data_r3,file="skmel475_r3.txt",col.names=F,row.names=F,quote=F,append=T)
write.table(t(head),file="skmel475_raw_data.txt",col.names=F,row.names=F,quote=F)

write.table(data_formatted,file="skmel475_raw_data.txt",col.names=F,row.names=F,quote=F,append=T)


#plot replicates and the correlation coeff. between replicates
cc12 <- cor(data_r1[,9],data_r2[,9])
pdf("replicates_12.pdf")
plot(data_r1[,9],data_r2[,9],xlab="replicate 1",ylab="replicate 2")
title(main=paste("skmel475 replicates \nCC=",cc12))

dev.off()

excl_list <- union(grep("6",data_r1[,3]),intersect(grep("DMSO",data_r1[,2]),grep("24",data_r1[,3])))
#print(excl_list)
data_r1_ex <-data_r1[-excl_list,]
data_r2_ex <-data_r2[-excl_list,]
pdf("replicates_13.pdf")
cc13 <- cor(data_r1_ex[,9],data_r3[,9])
#print(cc13)
plot(data_r1_ex[,9],data_r3[,9],xlab="replicate 1",ylab="replicate 3")
title(main=paste("skmel475 replicates \nCC=",cc13))

dev.off()

pdf("replicates_23.pdf")
cc23 <- cor(data_r2_ex[,9],data_r3[,9])
#print(cc13)
plot(data_r2_ex[,9],data_r3[,9],xlab="replicate 2",ylab="replicate 3")
title(main=paste("skmel475 replicates \nCC=",cc23))

dev.off()

 


rm(list=ls())  






