#plot proteomics data
library(Hmisc)
library(graphics)
nab=72
ntime=4
ncond=3
cond1="DMSO"
cond2="901"
cond3="4032"
cond4="combo"
data_r <- read.table("skmel_475_a_n_tp.txt",header=T)
abs <-read.table("ab_table.txt")

data_c1 <- data_r[grep(cond1,data_r[,2]),]
data_c2 <- data_r[grep(cond2,data_r[,2]),]
data_c3 <- data_r[grep(cond3,data_r[,2]),]
data_c4 <- data_r[grep(cond4,data_r[,2]),]
pdf("skmel475_panel4.pdf")
#order_data

par(mfrow=c(6,3),mar = rep(2, 4))
for(i in 55:72) {
	time = c(0,1,6,24)
#assume dmso(t=1)=dmso(t=0)	
abindex1=which(data_c1[,5]==abs[i,1])
dmso1_0=which((data_c1[,3]==1))

abindex2=which(data_c2[,5]==abs[i,1])
dmso2_0=which((data_c2[,3]==1))

abindex3=which(data_c3[,5]==abs[i,1])
dmso3_0=which((data_c3[,3]==1))

abindex4=which(data_c4[,5]==abs[i,1])
dmso4_0=which((data_c4[,3]==1))

readout1 = append(data_c1[intersect(dmso1_0,abindex1),9],data_c1[abindex1,9])     
err1=append(data_c1[intersect(dmso1_0,abindex1),12],data_c1[abindex1,12])
#print(as.matrix(data_c4[abindex1,]))

readout2 = append(data_c1[intersect(dmso1_0,abindex1),9],data_c2[abindex2,9])     
err2=append(data_c1[intersect(dmso1_0,abindex1),12],data_c2[abindex2,12])

readout3 = append(data_c1[intersect(dmso1_0,abindex1),9],data_c3[abindex3,9])     
err3=append(data_c1[intersect(dmso1_0,abindex1),12],data_c3[abindex3,12])

readout4 = append(data_c1[intersect(dmso1_0,abindex1),9],data_c4[abindex4,9])     
err4=append(data_c1[intersect(dmso1_0,abindex1),12],data_c4[abindex4,12])

opts = c("l")
tit=paste(abs[i,2],abs[i,3],sep="_")

#extension=".pdf"
#pdf(paste(toString(title, width=NULL),extension))
errbar(time,readout1,readout1+err1,readout1-err1,xlim=c(0,24),ylim=c(0,7),col=rgb(0,0,0),type=opts[1])
par(new=T)
errbar(time,readout2,readout2+err2,readout2-err2,xlim=c(0,24),ylim=c(0,7),col=rgb(1,0,0),ann=FALSE,type=opts[1])

par(new=T)
errbar(time,readout3,readout3+err3,readout3-err3,xlim=c(0,24),ylim=c(0,7),col=rgb(0,1,0),ann=FALSE,type=opts[1])

par(new=T)
errbar(time,readout4,readout4+err4,readout4-err4,xlim=c(0,24),ylim=c(0,7),col=rgb(0,0,1),,type=opts[1],xlab="time")
title(main=tit,cex.main=0.7)

#dev.off()
}
dev.off()
