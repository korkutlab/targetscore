nab=72
ncond=4
ntime=3
normfile="cellmix.txt"
totalprot <- read.table("total_prot_skmel475.txt",header=T) 

normal <- read.table(normfile)
data_r1 <- read.table("skmel475_r1.txt",header=T)
data_r2 <- read.table("skmel475_r2.txt",header=T)
data_r3 <- read.table("skmel475_r3.txt",header=T)

###normalize wrt total protein levels####
data_r1_tp <- as.matrix(data_r1)
data_r2_tp <- as.matrix(data_r2)
data_r3_tp <- as.matrix(data_r3)
totalprot_A <-totalprot[grep("A",totalprot[,2]),]
totalprot_B <-totalprot[grep("B",totalprot[,2]),]
totalprot_C <-totalprot[grep("C",totalprot[,2]),]
#print(totalprot_C[1,])
for(i in 1:nab){
	 for (j in 1:12){
	data_r1_tp[ncond*ntime*(i-1)+j,9]<-2*as.numeric(data_r1[ncond*ntime*(i-1)+j,9])/as.numeric(totalprot_A[j,12])	
	data_r2_tp[ncond*ntime*(i-1)+j,9]<-2*as.numeric(data_r2[ncond*ntime*(i-1)+j,9])/as.numeric(totalprot_B[j,12])	
	 }
		 for (j in 1:7){
		 	data_r3_tp[7*(i-1)+j,9]<-2*as.numeric(data_r3[7*(i-1)+j,9])/as.numeric(totalprot_C[j,12])
		 } 
}
print(as.numeric(data_r1_tp[,9])-as.numeric(data_r1[,9]))
############################################

data_m1 <- as.matrix(data_r1_tp)
data_m2 <- as.matrix(data_r2_tp)
data_m3 <- as.matrix(data_r3_tp)
#normalization data
norm1 <-as.matrix(normal[grep("32",normal[,7]),])

#average data
excl_list <- union(grep("6",data_r1[,3]),intersect(grep("DMSO",data_r1[,2]),grep("24",data_r1[,3])))
print(excl_list)
data_r1_inc <-(data_m1[-excl_list,])
data_r2_inc <-(data_m2[-excl_list,])
data_r1_ex <-(data_m1[excl_list,])
data_r2_ex <-(data_m2[excl_list,])
data_ave <- data_m1
a1a2<-data_r1_ex[,9]
a1a2a3<-data_r1_inc[,9]
a1m <- data_m1[,9]


data_ave[,4]<- "aver"
data_ave[excl_list,9] <- 0.5*(as.numeric(data_m1[excl_list,9])+as.numeric(data_m2[excl_list,9]))
print(length(excl_list))

#for (i in 1:length(excl_list)){
#	sd[excl_list]=sd((data_m1[excl_list,9],data_m2[excl_list,9]),na.rm=F)/sqrt(2)
#}
#for (i in 1:length(-excl_list)){
#	a1a2a3[i]=sd(c(data_r1_inc[i,9],data_r2_inc[i,9],data_m3[i,9]))/sqrt(3)
#}
#print(a1a2a3)
data_ave[-excl_list,9] <-(as.numeric(data_m1[-excl_list,9])+as.numeric(data_m2[-excl_list,9])+as.numeric(data_m3[,9]))/3

a1m[-excl_list] <- ((as.numeric(data_m1[-excl_list,9])-as.numeric(data_ave[-excl_list,9]))^2+(as.numeric(data_m2[-excl_list,9])-as.numeric(data_ave[-excl_list,9]))^2+(as.numeric(data_m3[,9])-as.numeric(data_ave[-excl_list,9]))^2)/3
a1m[excl_list] <- ((as.numeric(data_m1[excl_list,9])-as.numeric(data_ave[excl_list,9]))^2+(as.numeric(data_m2[excl_list,9])-as.numeric(data_ave[excl_list,9]))^2)/2
print(a1a2a3)
data_ave_c <-cbind(data_ave,a1m)
write.table(data_ave_c,file="skmel_475_a_tp.txt", row.names=F, quote=F)
#print(as.numeric(data_m1[-excl_list,9])+as.numeric(data_m3[,9]))
#data_ave_n <- data_ave 
data_ave_n_m <- (data_m1)
data_ave_n_m[,10]=0
sem<-(a1m)
#print(a1m)
#print(data_ave[-excl_list,11])

#print(data_ave[excl_list,9])
###################################################

for (i in 1:nab){
	for (j in 1:ntime){
	 for (k in 1:ncond){
#	 	print(k)
#        print(data_ave[(i-1)*ntime*ncond+(j-1)*ncond+k,9])
	 	data_ave_n_m[(i-1)*ntime*ncond+(j-1)*ncond+k,9]<-as.numeric(data_ave[(i-1)*ntime*ncond+(j-1)*ncond+k,9])/as.numeric(norm1[i,14])
#	    data_ave_n_m[(i-1)*ntime*ncond+(j-1)*ncond+k,11]<-as.numeric(data_ave[(i-1)*ntime*ncond+(j-1)*ncond+k,11])/as.numeric(norm1[i,14])	
        sem[(i-1)*ntime*ncond+(j-1)*ncond+k]=as.numeric(a1m[(i-1)*ntime*ncond+(j-1)*ncond+k])/as.numeric(norm1[i,14])	
#        print(norm1[i,14])
#        print(norm1[i,9])
         
	 }		
	}
}
#print(as.numeric(norm1[]))
data_ave_n_m_c <-cbind(data_ave_n_m,sem)

write.table(data_ave_n_m_c ,file="skmel_475_a_n_tp.txt", row.names=F, quote=F)
rm(list=ls())  


