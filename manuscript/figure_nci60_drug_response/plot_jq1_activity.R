library(ggrepel)
library(ggplot2)

a <- read.table("~/Downloads/dataset_a.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
b <- read.table("~/Downloads/dataset_b.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(a) <- paste0("V", 1:8)
colnames(b) <- paste0("V", 1:8)

n <- intersect(a[,1], b[,1])
ia <- which(a$Cell.Line %in% n)
ib <- which(b$Cell.Line %in% n)

x <- a[a[,1] %in% n, 3]
y <- b[b[,1] %in% n, 3]
cor.test(x, y)

plot(x, y)

c1 <- grepl("OV:", a[,1])
c2 <- grepl("OV:", b[,1])

head(a[c1, c(1,3)])
head(b[c2, c(1,3)])

plot(a[c1, 3], b[c2, 3])

#ggplot(dat, aes(x=xvar, y=yvar)) + geom_point() 

d <- merge(a[c1, c(1,3)], b[c2, c(1,3)], by.x="V1", by.y="V1")
colnames(d) <- c("line", "public", "private")

p <- ggplot(d)
p <- p + geom_point(aes(public, private))
p <- p + geom_text_repel(aes(public, private, label = line))
p <- p + theme_bw()
p

e <- merge(a[, c(1,3)], b[, c(1,3)], by.x="V1", by.y="V1")
colnames(e) <- c("line", "public", "private")

cor.test(e$public, e$private, method = "spearman")

p <- ggplot(e)
p <- p + geom_point(aes(public, private))
p <- p + geom_text_repel(aes(public, private, label = line))
p <- p + theme_bw()
p
