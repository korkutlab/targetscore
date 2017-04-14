#bt474 consensus
#consensus1 24h 3D resistance
#consensus2 48h 3D response
#consensus3 24h 3D resistance
#consensus448h 3D response
cellline <- "bt474"
drug1 <- "JQ1"
drug2 <- "iBET151"
drug3 <- "iBET726"
drug4 <- "iBET762"
time1<- "24h"
dimension <- "3D"
folder <- cellline
proteins <- read.table(file="protlist1.txt",header=F,row.names = 1)
list1 <- read.table(file=
         paste(cellline,"/TS_",time1,"_",dimension,"_",drug1,"_",cellline,".txt",sep=""),
         header=T,row.names = 1)
list2 <- read.table(file=
                      paste(cellline,"/TS_",time1,"_",dimension,"_",drug2,"_",cellline,".txt",sep=""),
                    header=T,row.names = 1)
list3 <- read.table(file=
                      paste(cellline,"/TS_",time1,"_",dimension,"_",drug3,"_",cellline,".txt",sep=""),
                    header=T,row.names = 1)
list4 <- read.table(file=
                      paste(cellline,"/TS_",time1,"_",dimension,"_",drug4,"_",cellline,".txt",sep=""),
                    header=T,row.names = 1)
score <- array(0,dim=c(nrow(list1),8))

rownames(score)<- proteins[,1]
for(i in 1:nrow(list1)){
  score[i,1] <- list1$FDRp[which(rownames(proteins)[i]==rownames(list1))]
  score[i,2] <- list2$FDRp[which(rownames(proteins)[i]==rownames(list2))]
  score[i,3] <- list3$FDRp[which(rownames(proteins)[i]==rownames(list3))]
  score[i,4] <- list4$FDRp[which(rownames(proteins)[i]==rownames(list4))]
  
  score[i,5] <- list1$change[which(rownames(proteins)[i]==rownames(list1))]
  score[i,6] <- list2$change[which(rownames(proteins)[i]==rownames(list2))]
  score[i,7] <- list3$change[which(rownames(proteins)[i]==rownames(list3))]
  score[i,8] <- list4$change[which(rownames(proteins)[i]==rownames(list4))]
}

for(i in 1:nrow(list1)){
  for(j in 1:4){
    if(score[i,j] > 0.25 || any(is.na(score[i,j]))){score[i,j]<-0}      
      
  }
}
for(i in 1:nrow(list1)){
  for(j in 1:4){
    if(score[i,j] != 0 &(score[i,j+4]>0)){score[i,j]<- 1}      
    
  }
}
for(i in 1:nrow(list1)){
  for(j in 1:4){
    if(score[i,j] != 0 &(score[i,j+4]<0)){score[i,j]<- -1}      
    
  }
}
scoreT <- array(0,dim=c(nrow(list1)))
scoreT <- score[,1]+score[,2]+score[,3]+score[,4]
write.table(score[,1:4],file=paste("scoreTS_",time1,"_",dimension,"_",cellline,".txt",sep=""),quote=F,row.names=T)

resistance <- array(0,dim=c(length(proteins[which(scoreT >= 2),1]),4))
colnames(resistance) <- c("ID","include","change","FDRp")
resistance[,1]<- as.vector(proteins[which(scoreT >= 2),1])
resistance[,2]<-0
resistance[,3]<-5
resistance[,4]<-0.25
write.table(resistance,file=paste("resistConsTS_",time1,"_",dimension,"_",cellline,".txt",sep=""),quote=F,row.names=F)

#response
response <- array(0,dim=c(length(proteins[which(scoreT <= -2),1]),4))
colnames(response) <- c("ID","include","change","FDRp")
response[,1]<- as.vector(proteins[which(scoreT <= -2),1])
response[,2]<-0
response[,3]<--5
response[,4]<-0.25
write.table(response,file=paste("responseConsTS_",time1,"_",dimension,"_",cellline,".txt",sep=""),quote=F,row.names=F)
rm(list = ls())
