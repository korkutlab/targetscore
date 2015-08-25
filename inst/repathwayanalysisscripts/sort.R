nab=72
nrep=2
ncond=3
ntime=4
ncell_line=2

data <- read.table("zepto_data_skmel118_475.txt",header=T)
print(data)
#dataqual <- as.character(data[,6])

datastr <- as.character(data[,8])
#dmel <- grep("Melanoma",datastr)
dnorm <-grep("Mix",datastr)
write.table(data[dnorm,], file="cellmix.txt", col.names=F, row.names=F, sep = "	",quote=FALSE)

#write.table(data[dmel,], file="skmel475.txt", col.names=F, row.names=F, sep = "	",quote=FALSE)

dskmel475 <- grep("Mel-475",datastr)
dskmel118 <- grep("Mel-118",datastr)
write.table(data[dskmel475,], file="skmel475.txt", col.names=F, row.names=F, sep = "	",quote=FALSE)
write.table(data[dskmel118,], file="skmel118.txt", col.names=F, row.names=F, sep = "	",quote=FALSE)

dataskmel475=data[dskmel475,]
dataskmel118=data[dskmel118,]






