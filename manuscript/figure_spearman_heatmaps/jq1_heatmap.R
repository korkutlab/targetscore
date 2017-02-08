library(ggplot2)
library(reshape2)
library(gridExtra)

rm(list = ls())

# Based on https://rpubs.com/danmirman/plotting_factor_analysis
workDir <- "figure_spearman_heatmaps"
basePath <- file.path("figure_spearman_heatmaps", "spearmanData")
cellLines <- c("cov318", "igrov1", "ovcar3", "ovcar4")

grobs <- list()

for(cellLine in cellLines) {
  #cellLine <- "cov318"

  heatmapDat <- read.table(paste0(basePath, cellLine, "_jq1_ave.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
  sprDat <- read.table(paste0(basePath, cellLine, "_spr_sorted.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

  # DEBUG
  #sprDat <- data_frame(name=colnames(heatmapDat[,2:51]), value=rnorm(50))

  dropAbnames <- sprDat[sprDat$adjPval > 0.05, "abName"]
  sprDat <- sprDat[sprDat$adjPval <= 0.05, ]

  # Force ggplot order. Even though decreasing=FALSE it is actually decreasing
  sprDat <- sprDat[order(sprDat[, "corr"], decreasing=FALSE), ]
  sortedAbName <- sprDat$abName
  sprDat$abName <- factor(sprDat$abName, levels=sprDat$abName)

  dropCol <- c("GORDER", dropAbnames)

  # NOTE: One extra column that has the antibody names
  heatmapDat <- heatmapDat[, !(names(heatmapDat) %in% dropCol)]

  meltedDat <- melt(heatmapDat)

  meltedDat$variable <- factor(meltedDat$variable, levels=sortedAbName)

  #breaks <- c(-2, -1, 0, 1, 2)

  # Heatmap and barplot
  # From: https://rpubs.com/danmirman/plotting_factor_analysis
  p1 <- ggplot(data=meltedDat, aes_string(x="Sample_description", y="variable", fill="value")) +
    geom_tile() +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint=0, limit=c(-2, 2), space = "Lab",
                         na.value = 'green',
                         name="Pearson\nCorrelation") +
    theme_minimal() +
    theme(
      plot.margin=unit(c(2,2,10.2,2), "mm"),
      plot.background=element_blank(),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      axis.text.x=element_blank(),
      #axis.text.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position = "none"
    ) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ggtitle(toupper(cellLine)) + theme(plot.title=element_text(hjust=0.5))
    #guides(fill=FALSE) #omit unnecessary gradient legend

  p2 <- ggplot(sprDat, aes_string(x="abName", y="corr")) +
    geom_bar(stat="identity", fill="gray70") +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.margin=unit(c(7.5,3,3,3), "mm"),
      plot.background=element_blank(),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      #axis.text.x=element_text(angle = -90, hjust=1),
      #axis.text.x=element_blank(),
      #axis.text.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks=seq(-1, 1, 1))
    #+
    #scale_x_discrete(expand=c(0,0)) +
    #scale_y_discrete(expand=c(0,0))

  # Hide the barplot labels
  # NOTE: You can overwrite individual theme entries
  p2 <- p2 + theme(axis.text.y=element_blank())
  #p2

  t1 <- paste0(cellLine, "_heatmap")
  t2 <- paste0(cellLine, "_bar")

  #gA <- ggplot_gtable(ggplot_build(p1))
  #gB <- ggplot_gtable(ggplot_build(p2))

  #grobs[[t1]] <- gA
  #grobs[[t2]] <- gB

  grobs[[t1]] <- ggplotGrob(p1)
  grobs[[t2]] <- ggplotGrob(p2)
}

# NOTE: The legend was brought in from a separate plot.

# Plot the various elements (e.g. heatmaps/bar plots) on a single plot
base <- c(1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8)
#colorBar <- rep(9, length(base))
#base <- c(1,2,3,4,5,6,7,8)

lay <- rbind(base, base, base)

#gA <- ggplot_gtable(ggplot_build(p1))
#gB <- ggplot_gtable(ggplot_build(p2))

#gA <- grobs[[t1]]
#gB <- grobs[[t2]]
#grid.arrange(gA, gB, gC, gD, nrow=1, layout_matrix=lay)
#grid.arrange(gA, gB, gA, gB, gA, gB, gA, gB, nrow=1, layout_matrix=lay)

#??? NO IDEA
#grid.arrange(grobs, nrow=1, layout_matrix=lay)

#grob <- arrangeGrob(gA, gB, gA, gB, gA, gB, gA, gB, nrow=1, layout_matrix=lay)
grob <- arrangeGrob(grobs=grobs, nrow=1, layout_matrix=lay)
ggsave(file=file.path(workDir, "jq1_spearman_heatmap.pdf"), grob, width=15, height=7.64)
