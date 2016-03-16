#data <- read.table("skmel475.txt")


#ab <- data[grep("_DMSO_1_A",data[,8]),9]

#write.table(ab,file="Ab_index.txt",col.names=F,row.names=F,quote=F)
ab <- read.table("ab_index.txt")
ab_o <- read.table("ab_index_o.txt")

abx <- as.matrix(ab)
print(abx)
for (i in 1:72 ){
	 for (k in 1:72){	
	 	if(as.character(ab[i,1])==as.character(ab_o[k,1])){
#	 		print(ab[i])
	abx[i,]=as.matrix(ab_o[k,])

	}	
	
	}

	
}
write.table(abx, file="match.txt",col.names=F,row.names=F,quote=F)

