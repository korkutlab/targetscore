nab=72
ncond=4
ntime=3
normfile="cellmix.txt"
normal <- read.table(normfile)
data_r1 <- read.table("skmel475_r1.txt",header=T)
data_r2 <- read.table("skmel475_r2.txt",header=T)
data_r3 <- read.table("skmel475_r3.txt",header=T)
data_m1 <- as.matrix(data_r1)
data_m2 <- as.matrix(data_r2)
data_m3 <- as.matrix(data_r3)

#average data
excl_list <- union(grep("6",data_r1[,3]),intersect(grep("DMSO",data_r1[,2]),grep("24",data_r1[,3])))
#print(excl_list)
data_r1_ex <-(data_m1[-excl_list,])
data_r2_ex <-(data_m2[-excl_list,])
data_r1_inc <-(data_m1[excl_list,])
data_r2_inc <-(data_m2[excl_list,])
data_ave <- data_m1
data_ave[,4]<- "aver"
data_ave[excl_list,9] <- 0.5*(as.numeric(data_m1[excl_list,9])+as.numeric(data_m2[excl_list,9]))
data_ave[-excl_list,9] <-(as.numeric(data_m1[-excl_list,9])+as.numeric(data_m2[-excl_list,9])+as.numeric(data_m3[,9]))/3
write.table(data_ave,file="skmel_475_a.txt", row.names=F, quote=F)
#print(as.numeric(data_m1[-excl_list,9])+as.numeric(data_m3[,9]))
#data_ave_n <- data_ave 
data_ave_n_m <- (data_ave)
#normalization data
norm1 <-as.matrix(normal[grep("64",normal[,7]),])
#print(norm1)



############################################
for (i in 1:nab){
	for (j in 1:ntime){
	 for (k in 1:ncond){
#	 	print(k)
#        print(data_ave[(i-1)*ntime*ncond+(j-1)*ncond+k,9])
	 	data_ave_n_m[(i-1)*ntime*ncond+(j-1)*ncond+k,9]<-as.numeric(data_ave[(i-1)*ntime*ncond+(j-1)*ncond+k,9])/as.numeric(norm1[i,14])
	 }		
	}
}

write.table(data_ave_n_m,file="skmel_475_a_n.txt", row.names=F, quote=F)



rm(list=ls())  


