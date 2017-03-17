library(heatmap3)
#TS_cov318_24h_3D_JQ1 <- read.table("ts_cov318/cov318_JQ1_3D_24h_TS.txt",header=F,row.names=1)
#TS_ovcar3_24h_3D_JQ1 <- read.table("ts_ovcar3/ovcar3_JQ1_3D_24h_TS.txt",header=F,row.names=1)
#TS_ovcar4_24h_3D_JQ1<- read.table("ts_ovcar4/ovcar4_JQ1_3D_24h_TS.txt",header=F,row.names=1)
TS_skov3_24h_3D_JQ1<- read.table("ts_skov3/skov3_JQ1_3D_24h_TS.txt",header=F,row.names=1)
TS_bt474_24h_3D_JQ1 <- read.table("ts_bt474/bt474_JQ1_3D_24h_TS.txt",header=F,row.names=1)
TS_hcc1954_24h_3D_JQ1 <- read.table("ts_hcc1954/hcc1954_JQ1_3D_24h_TS.txt",header=F,row.names=1)
TS_mdamb468_24h_3D_JQ1 <- read.table("ts_mdamb468/mdamb468_JQ1_3D_24h_TS.txt",header=F,row.names=1)
TS_skbr3_24h_3D_JQ1 <- read.table("ts_skbr3/skbr3_JQ1_3D_24h_TS.txt",header=F,row.names=1)

#q_cov318_24h_3D_JQ1 <- read.table("ts_cov318/cov318_JQ1_3D_24h_q.txt")
#q_ovcar3_24h_3D_JQ1 <- read.table("ts_ovcar3/ovcar3_JQ1_3D_24h_q.txt")
#q_ovcar4_24h_3D_JQ1<- read.table("ts_ovcar4/ovcar4_JQ1_3D_24h_q.txt")
q_skov3_24h_3D_JQ1<- read.table("ts_skov3/skov3_JQ1_3D_24h_q.txt")
q_bt474_24h_3D_JQ1 <- read.table("ts_bt474/bt474_JQ1_3D_24h_q.txt")
q_hcc1954_24h_3D_JQ1 <- read.table("ts_hcc1954/hcc1954_JQ1_3D_24h_q.txt")
q_mdamb468_24h_3D_JQ1 <- read.table("ts_mdamb468/mdamb468_JQ1_3D_24h_q.txt")
q_skbr3_24h_3D_JQ1 <- read.table("ts_skbr3/skbr3_JQ1_3D_24h_q.txt")
TS_bt474_24h_3D_JQ1 <- as.data.frame(TS_bt474_24h_3D_JQ1[-36,])
colnames(TS_bt474_24h_3D_JQ1)<-colnames(TS_mdamb468_24h_3D_JQ1)
rownames(TS_bt474_24h_3D_JQ1)<-rownames(TS_mdamb468_24h_3D_JQ1)
TS_24h_3D_JQ1 <- cbind(TS_skov3_24h_3D_JQ1,
                         TS_bt474_24h_3D_JQ1[-36,],
                         TS_hcc1954_24h_3D_JQ1,
                         TS_mdamb468_24h_3D_JQ1,
                         TS_skbr3_24h_3D_JQ1[-36,])
colnames(TS_24h_3D_JQ1) <- c("SKOV3","BT474","HCC1954","MDAMB468","SKBR3")

med_TS_24h_3D_JQ1 <- apply(TS_24h_3D_JQ1, 1, median)
TS_24h_3D_JQ1$median <-med_TS_24h_3D_JQ1
TS_24h_3D_JQ1_Srt <- TS_bt474_24h_3D_JQ1[order(TS_bt474_24h_3D_JQ1, decreasing = T),]
write.table(TS_24h_3D_JQ1_Srt,file="TS_24h_3D_JQ1.txt",quote=F)
#TS_bt474_24h_3D_JQ1 <- TS_bt474_24h_3D_JQ1[-36,]
pdf("3D.pdf")
heatmap3(TS_24h_3D_JQ1_Srt,                     
          Rowv=NA, Colv=NA,          
          cexRow=0.5,cexCol=0.75)     
dev.off()
rm(list = ls())
