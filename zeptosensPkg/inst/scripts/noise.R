##use median of each condition instead of all replicates
data <- read.table("skmel475_raw_data.txt", header=T)
#print(data[,5])
ab_index <- read.table("ab_table.txt")
secab <-grep("SECAB",data[,6])
nab=72
sd_secab <- sd(data[secab,9])
ave_secab<- mean(data[secab,9])
max_secab<- max(data[secab,9])
min_secab<- min(data[secab,9])
#print(sd_secab)
#print(ave_secab)
#print(max_secab)
#print(min_secab)
for (i in 1:nab){
#	print(data[which(data[,5]==ab_index[i,1]),5])
	maxval_ab=max(data[which(data[,5]==ab_index[i,1]),9])
#	print(maxval_ab)
	if(maxval_ab < (as.numeric(ave_secab)+3*as.numeric(sd_secab ))){
		print(ab_index[i,1])
#print(as.matrix(data[grep(as.character(ab_index[i,1]),(data[,5])),5]))
		
		
	}
}

rm(list=ls())  

