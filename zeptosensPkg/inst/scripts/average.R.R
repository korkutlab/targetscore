data1 <- read.table("skmel475_r1.txt")
data2 <- read.table("skmel475_r2.txt")
data3 <- read.table("skmel475_r3.txt")

#average data
excl_list <- union(grep("6",data_r1[,3]),intersect(grep("DMSO",data_r1[,2]),grep("24",data_r1[,3])))
#print(excl_list)
data_r1_ex <-data_r1[-excl_list,]
data_r2_ex <-data_r2[-excl_list,]

data_ave[excl_list,9] <- 0.5*(data_r1[excl_list,9]+data_r2[excl_list,9])
data_ave[-excl_list,9] <-(data_r1[-excl_list,9]+data_r2[-excl_list,9]+data_r3[,9])/3

print(data_ave)

