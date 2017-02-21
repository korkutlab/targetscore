options("digits"=3)
nrep=3
ncond=3
ntime=4
nab=72
ctr=ncond*ntime*nrep
ct=ncond*ntime
data <- read.table("skmel475.txt")
annot <- read.table("ab_table.txt")
normal <- read.table("cellmix.txt")
data_r1 <- data[grep("_A",data[,8]),]
data_r2 <- data[grep("_B",data[,8]),]
data_r3 <- data[grep("_C",data[,8]),]
#print(data_1)
write.table(data_r1,file="skmel475_r1.txt",col.names=F,row.names=F,quote=F)
write.table(data_r2,file="skmel475_r2.txt",col.names=F,row.names=F,quote=F)
write.table(data_r3,file="skmel475_r3.txt",col.names=F,row.names=F,quote=F)
cc1<-cor(data_r1[,14],data_r2[,14])
print(cc1)

plot(data_r1[,14],data_r2[,14])
excl_list <- union(grep("_6_",data_r1[,8]),grep("DMSO_24",data_r1[,8]))
#print(excl_list)
data_r1_ex <-data_r1[-excl_list,]
data_r2_ex <-data_r2[-excl_list,]
cc2<-cor(data_r1_ex[,14],data_r3[,14])
print(cc2)
#plot(data_r1_ex[,14],data_r3[,14])
cc3<-cor(data_r2_ex[,14],data_r3[,14])
print(cc3)
#plot(data_r2_ex[,14],data_r3[,14])


#(for(i in 1:nab)((i-1)*ncond*ntime+2)),]
#cc2 <- cor(data_r1[(grep(,data_r3[])),14],data_r3[,14])
#print(cc1)
#bb <- intersect(grep("_A",data[,8]) , grep("Akt",data[,9]))
#for (i in 1:nab){
#	for (k in 1:ncond){
#		for (j in 1:ntime)
#		fs[i]
#		
#	}
#}

#print(bb)

