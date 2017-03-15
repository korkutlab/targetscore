library(ggplot2)

source("https://gist.githubusercontent.com/cannin/db5174a74349e601fbcd662f1fa2951f/raw/saveGgplotPlot.R")

#tmpDat has 2 columns FC, negLogPadj

p <- ggplot() +
  geom_point(data=tmpDat, aes(x=FC, y=negLogPadj, color=color, text=toolTip), alpha=0.4, size=2) +
  theme_bw() +
  xlab("log2 (Tumor/Normal)") + ylab("-log10 P-Value") +
  scale_color_manual(name="", values=c("red", "black"))

saveGgplotPlot(p)
