# Basic plot
# plot(ovcar3_TS_all$X4, -ovcar3_TS_all$X3)
# cor(ovcar3_TS_all$X4, -ovcar3_TS_all$X3, use = "pairwise.complete.obs", method = "spearman")
# abline(lm(X3 ~ X4, ovcar3_TS_all))

library(RColorBrewer)
library(ggplot2)

ptCol <- brewer.pal(8, "Set1")[2]
ptCol <- "black"

ovcar3_TS_all$negCorr <- -ovcar3_TS_all$X3
p <- ggplot(ovcar3_TS_all, aes(negCorr, X4))
p <- p + geom_point(colour = ptCol) + ggtitle("Target Score") + xlab("Dose Dependent Response") + ylab("Target Score") + geom_hline(yintercept=0) + geom_vline(xintercept=0) + theme_bw()
p


