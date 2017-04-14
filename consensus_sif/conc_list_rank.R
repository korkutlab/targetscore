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
score <- array(0,dim=c(nrow(list1),12))

rownames(score)<- rownames(proteins)
for(i in 1:nrow(list1)){
  score[i,1] <- list1$change[which(rownames(proteins)[i]==rownames(list1))]
  score[i,2] <- list2$change[which(rownames(proteins)[i]==rownames(list2))]
  score[i,3] <- list3$change[which(rownames(proteins)[i]==rownames(list3))]
  score[i,4] <- list4$change[which(rownames(proteins)[i]==rownames(list4))]
  
}

score[,5] <- match(score[,1],sort(score[,1])) #list1$change[which(rownames(proteins)[i]==rownames(list1))]
score[,6] <- match(score[,2],sort(score[,2]))#list2$change[which(rownames(proteins)[i]==rownames(list2))]
score[,7] <- match(score[,3],sort(score[,3]))#list3$change[which(rownames(proteins)[i]==rownames(list3))]
score[,8] <- match(score[,4],sort(score[,4]))#list4$change[which(rownames(proteins)[i]==rownames(list4))]
for(i in 1:nrow(list1)){
  if(score[i,5]<11){score[i,9]<-  -1} 
  if(score[i,6]<11){score[i,10]<- -1} 
  if(score[i,7]<11){score[i,11]<- -1} 
  if(score[i,8]<11){score[i,12]<- -1} 

  if(score[i,5]>207){score[i,9]<-  1} 
  if(score[i,6]>207){score[i,10]<- 1} 
  if(score[i,7]>207){score[i,11]<- 1} 
  if(score[i,8]>207){score[i,12]<- 1}   
}



scoreT <- array(0,dim=c(nrow(list1)))
scoreT <- score[,9]+score[,10]+score[,11]+score[,12]
write.table(score[,9:12],file=paste("scoreTS_",time1,"_",dimension,"_",cellline,".txt",sep=""),quote=F,row.names=T)

resistance <- array(0,dim=c(length(proteins[which(scoreT >= 3),1]),4))
colnames(resistance) <- c("ID1","include","change","FDRp")
resistance[,1]<- as.vector(rownames(proteins)[which(scoreT >= 3)])
resistance[,2]<-0
resistance[,3]<-5
resistance[,4]<-0.25
write.table(resistance,file=paste("resistConsTS_",time1,"_",dimension,"_",cellline,".txt",sep=""),quote=F,row.names=F,sep="\t")

#response
response <- array(0,dim=c(length(proteins[which(scoreT <= -3),1]),4))
colnames(response) <- c("ID1","include","change","FDRp")
response[,1]<- as.vector(rownames(proteins)[which(scoreT <= -3)])
response[,2]<-0
response[,3]<--5
response[,4]<-0.25
write.table(response,file=paste("responseConsTS_",time1,"_",dimension,"_",cellline,".txt",sep=""),quote=F,row.names=F,sep="\t")
rm(list = ls())
