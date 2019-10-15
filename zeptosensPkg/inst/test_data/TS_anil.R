library(plyr)
library(zeptosensPkg)
library(data.table)
library(RColorBrewer)
library(gplots)
library(matrixStats)
library(ggplot2)
library(glasso)
library(metap)
library(ggrepel)
library(pheatmap)

set.seed(1)

# FIXME fs.csv or fs.txt

# SET PARAMETERS ----
# data_dir <- "../data(ts)/ts_anil"
# protein_list_dir <- "../data(ts)/Protein_Name_List"
data_dir <- "inst/test_data"
resource_dir <- "inst/targetscoreData"
output_dir <- "inst/test_data/output"
max_dist <- 1 # changing this value requires additional work to compute product(wk). This is not a priority

# READ ANTIBODY FILE ----
antibody_map_file <- read.table(file.path(resource_dir, "antibodyMapFile.txt"),
  sep = "\t",
  header = TRUE, stringsAsFactors = FALSE
)

# READ IN DATA ----
sample1 <- "BT474"

# Read proteomic response for cellline1
proteomic_responses <- read.csv(file.path(data_dir, "BT474.csv"), row.names = 1) # n_prot=304

n_prot <- 304
proteomic_responses <- proteomic_responses
max_dist <- 1
antibody_map_file <- antibody_map_file
dist_file <- NULL
verbose <- FALSE

# Extract network for Ovarian Cancer from SignedPC
network <- zeptosensPkg::predict_bio_network(
  n_prot = dim(proteomic_responses)[2],
  proteomic_responses = proteomic_responses,
  max_dist = 1,
  antibody_map_file = antibody_map_file
)

wk <- network$wk # sum(wk!=0)=834
wks <- network$wks
dist_ind <- network$dist_ind
inter <- network$inter

# FIXME: NECESSARY? ADDED AL
write.csv(wk, file.path(output_dir, "wk.csv"))
write.csv(dist_ind, file.path(output_dir, "wk.csv"))

# adjust for fs.txt as the functional node
fs <- read.csv(file.path(resource_dir, "fs.csv"), header = TRUE, stringsAsFactors = FALSE)

# BT474 ----
# Calculate Target Score
n_prot <- ncol(proteomic_responses)
n_cond <- nrow(proteomic_responses)
target_score_output_file <- file.path(output_dir, "ts_bt474.txt")
matrix_wk_output_file <- file.path(output_dir, "wk_BT474.txt")
signedmatrix_wk_output_file <- file.path(output_dir, "wks_BT474.txt")

target_score_q_value_file <- file.path(output_dir, paste0(sample1, "_q.txt"))
target_score_dose_file <- file.path(output_dir, paste0(sample1, "_ts_d.txt"))
target_score_p_value_file <- file.path(output_dir, paste0(sample1, "_p.txt"))

ts <- array(0, dim = c(n_cond, n_prot))
ts_p <- array(0, dim = c(n_cond, n_prot))
ts_q <- array(0, dim = c(n_cond, n_prot))

n_perm <- 25

for (i in 1:n_cond) {
  results <- zeptosensPkg::get_target_score(
    wk = wk,
    wks = wks,
    dist_ind = dist_ind,
    inter = inter,
    n_dose = 1,
    n_prot = n_prot,
    proteomic_responses = proteomic_responses[i, ],
    max_dist = max_dist,
    n_perm = n_perm,
    cell_line = sample1,
    verbose = FALSE,
    fs_file = file.path(resource_dir, "fs.txt"),
    target_score_output_file = target_score_output_file,
    matrix_wk_output_file = matrix_wk_output_file,
    target_score_q_value_file = target_score_q_value_file,
    target_score_dose_file = target_score_dose_file,
    target_score_p_value_file = target_score_p_value_file
  )

  ts[i, ] <- results$ts
  ts_p[i, ] <- results$pts
  ts_q[i, ] <- results$q
}

# FIXME: WHY THIS?
colnames(ts) <- colnames(proteomic_responses)
ts <- data.frame(rownames(proteomic_responses), ts)
write.csv(ts, file = file.path(output_dir, "ts_bt474.csv"), row.names = FALSE)

colnames(ts_p) <- colnames(proteomic_responses)
ts_p <- data.frame(rownames(proteomic_responses), ts_p)
write.csv(ts_p, file = file.path(output_dir, "ts_bt474_pvalue.csv"), row.names = FALSE)

colnames(ts_q) <- colnames(proteomic_responses)
ts_q <- data.frame(rownames(proteomic_responses), ts_q)
write.csv(ts_q, file = file.path(output_dir, "ts_bt474_qvalue.csv"), row.names = FALSE)

# FIXME dist_ind.txt, randts.txt

# PLOT ----
# Heatmap
ts_bt474 <- read.csv(file.path(output_dir, "ts_bt474.csv"), row.names = 1)
data <- ts_bt474

bk <- c(seq(-4, -0.1, by = 0.01), seq(0, 5.5, by = 0.01))

# FIXME
filename <- file.path(output_dir, "heatmap_BT474_pdf")
pdf(filename)
pheatmap(data,
  scale = "none",
  color = c(
    colorRampPalette(colors = c("navy", "white"))(length(seq(-4, -0.1, by = 0.01))),
    colorRampPalette(colors = c("white", "firebrick3"))(length(seq(0, 5.5, by = 0.01)))
  ),
  legend_breaks = seq(-4, 5.5, 2), cellwidth = 1, cellheight = 1, fontsize = 1, fontsize_row = 1,
  breaks = bk
)
dev.off()

# Volcano Plot
## BT474
data <- read.csv(file.path(output_dir, "ts_bt474.csv"), row.names = 1)
data_p <- read.csv(file.path(output_dir, "ts_bt474_qvalue.csv"), row.names = 1)
data <- as.matrix(data)
data_p <- as.matrix(data_p)

for (i in 1:nrow(data)) {
  get_volcano_plot(ts = data[i, ], q_value = data_p[i, ], filename = rownames(data)[i], path = output_dir)
}

# Subnetwork for Top 30 Proteins ----
##  BT474
data <- read.csv(file.path(output_dir, "ts_bt474.csv"), row.names = 1)

# mcl1
data <- t(data[1, ])
data <- data.frame(data = data, prot = rownames(data))
data <- data[order(data$BT474_MCL1, decreasing = TRUE), ]
data_30 <- data.frame(prot = rownames(data)[1:30], ts = data$BT474_MCL1[1:30])

# WHY wk.csv read
wk <- read.csv(file = file.path(output_dir, "wk.csv"), row.names = 1)
index <- which(colnames(wk) %in% data_30$prot)
subnet <- wk[index, index]
edgelist_subnet <- zeptosensPkg::create_sif_from_matrix(t_net = subnet, genelist = colnames(subnet))
write.table(edgelist_subnet,
  file.path(output_dir, "signed_pc_BT474_MCL1_TOP30_positive_subnet.txt"),
  quote = FALSE,
  row.names = FALSE
)

data <- data[order(data$BT474_MCL1, decreasing = FALSE), ]
data_30 <- data.frame(prot = rownames(data)[1:30], ts = data$BT474_MCL1[1:30])

wk <- read.csv(file = file.path(output_dir, "wk.csv"), row.names = 1)
index <- which(colnames(wk) %in% data_30$prot)
subnet <- wk[index, index]
edgelist_subnet <- zeptosensPkg::create_sif_from_matrix(t_net = subnet, genelist = colnames(subnet))
write.table(edgelist_subnet,
  file.path(output_dir, "signed_pc_BT474_MCL1_TOP30_negative_subnet.txt"),
  quote = FALSE,
  row.names = FALSE
)

# COMB
data <- read.csv(file.path(output_dir, "ts_bt474.csv"), row.names = 1)
data <- t(data[2, ])
data <- data.frame(data = data, prot = rownames(data))
data <- data[order(data$BT474_comb, decreasing = TRUE), ]
data_30 <- data.frame(prot = rownames(data)[1:30], ts = data$BT474_comb[1:30])

wk <- read.csv(file = file.path(output_dir, "wk.csv"), row.names = 1)
index <- which(colnames(wk) %in% data_30$prot)
subnet <- wk[index, index]
edgelist_subnet <- zeptosensPkg::create_sif_from_matrix(t_net = subnet, genelist = colnames(subnet))
write.table(edgelist_subnet,
  file.path(output_dir, "signed_pc_BT474_COMB_TOP30_positive_subnet.txt"),
  quote = FALSE,
  row.names = FALSE
)

data <- data[order(data$BT474_comb, decreasing = FALSE), ]
data_30 <- data.frame(prot = rownames(data)[1:30], ts = data$BT474_comb[1:30])

wk <- read.csv(file = file.path(output_dir, "wk.csv"), row.names = 1)
index <- which(colnames(wk) %in% data_30$prot)
subnet <- wk[index, index]
edgelist_subnet <- zeptosensPkg::create_sif_from_matrix(t_net = subnet, genelist = colnames(subnet))
write.table(edgelist_subnet,
  file.path(output_dir, "signed_pc_BT474_COMB_TOP30_negative_subnet.txt"),
  quote = FALSE,
  row.names = FALSE
)
