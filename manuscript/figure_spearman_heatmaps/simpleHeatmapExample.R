library(ggplot2)
library(reshape2)
library(readr)

mydata <- mtcars[, c(1,3,4,5,6,7)]
cormat <- round(cor(mydata),2)

melted_cormat <- melt(cormat, na.rm=TRUE)

# Heatmap
p <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    legend.position = "bottom") +
  coord_fixed()
p 
ggsave("plot.pdf", p)

