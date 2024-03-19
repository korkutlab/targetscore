library(targetscore)

source("R/cv_glasso.R")

set.seed(123)

x <- matrix(runif(100), 5, 20)

cv_glasso(x)

verbose <- TRUE
max_dist <- 1

# NEW DATA: https://mdandersonorg-my.sharepoint.com/personal/hwang29_mdanderson_org/_layouts/15/onedrive.aspx?ga=1&id=%2Fpersonal%2Fhwang29%5Fmdanderson%5Forg%2FDocuments%2FBoxMigration%2FTS%2FTS%5FBC 

# outersect also known as symmetric difference 
outersect <- function(a,b) { setdiff(union(a,b), intersect(a,b)) }

data_dir <- system.file(file.path("extdata", "tcga_rppa"), package = "targetscore")
# TCGA-BRCA-L4.csv in test_data too small for predict_bio_network(); no interactions end up being found for parts of the code
#rppa_bc <- read.csv(file=file.path(data_dir, "TCGA-BRCA-L4.csv"), row.names=1) 
tmp <- read.csv(file=file.path(data_dir, "TCGA-BRCA-L4_v4.2.csv"), row.names=1)
rppa_bc_1 <- tmp[, 4:ncol(tmp)]

# Read drug perturbation data for BT474
data_dir <- system.file("test_data", package = "targetscore")
proteomic_responses <- read.csv(file.path(data_dir, "BT474.csv"), row.names = 1) # n_prot=61
tmp <- read.csv(file.path(data_dir, "cross_validation", "RPPA_HCC1954.csv"), row.names = 1)
rownames(tmp) <- paste(tmp$cellline, 
                       tmp$drug, 
                       tmp$drug_dose, 
                       tmp$time, 
                       tmp$Sample.description, 
                       sep="_")
proteomic_responses <- tmp[, 7:ncol(tmp)]

resource_dir <- system.file("target_score_data", package = "targetscore")
#mab_to_genes <- read.csv(file.path(resource_dir, "antibody_map_08092019.csv"), stringsAsFactors = FALSE)
mab_to_genes <- read.csv(file.path(resource_dir, "antibody_map_12162021.csv"), stringsAsFactors = FALSE)

network_ref <- targetscore::predict_bio_network(
  n_prot = ncol(proteomic_responses),
  proteomic_responses = proteomic_responses,
  max_dist = max_dist,
  mab_to_genes = mab_to_genes,
  verbose = verbose
)

# dephos_gene <-  match_genes_to_edgelist(genes1 = tmp_genes_ac, genes2 = tmp_genes_ai, annot_edgelist = dephosp, antibody_vec = colnames(proteomic_responses), use_annot = FALSE, verbose = verbose)
# 
# debugDat <- list(genes1=tmp_genes_ac, genes2=tmp_genes_ai, annot_edgelist=dephosp, antibody_vec = colnames(proteomic_responses), use_annot = FALSE, verbose = verbose); saveRDS(debugDat, "debug_dat.rds")
# 
# a <- readRDS("debug_dat.rds")
# genes1 <- a$genes1
# genes2 <- a$genes2
# annot_edgelist <- a$annot_edgelist
# antibody_vec <- a$antibody_vec
# use_annot <- a$use_annot
# verbose <- TRUE
# b <-  match_genes_to_edgelist(genes1 = genes1, genes2 = genes2, annot_edgelist = annot_edgelist, antibody_vec = antibody_vec, use_annot = use_annot, verbose = verbose)

# ERROR if using test_data rppa_bc not if using the full table
n_prot <- ncol(rppa_bc_1)
proteomic_responses <- rppa_bc_1
max_dist <- max_dist
mab_to_genes <- mab_to_genes

network_ref <- targetscore::predict_bio_network(
  n_prot = n_prot,
  proteomic_responses = proteomic_responses,
  max_dist = max_dist,
  mab_to_genes = mab_to_genes,
  verbose = FALSE
)
setdiff(colnames(proteomic_responses), mab_to_genes$AntibodyLabel)

# dist_file <- system.file("target_score_data", "distances.txt", package = "targetscore")
# tmp_dist <- read.table(dist_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 
# dist <- tmp_dist[which(tmp_dist[, 3] <= 1), ]
# idx_ab_map <- which(mab_to_genes[, 1] %in% colnames(proteomic_responses))
# 
# mab_genes <- mab_to_genes[idx_ab_map, 4]
# names(mab_genes) <- mab_to_genes[idx_ab_map, 1] 
# 
# dist_list <- match_genes_to_edgelist(
# genes1 = mab_genes,
# genes2 = NULL,
# annot_edgelist = dist,
# antibody_vec = colnames(proteomic_responses),
# use_annot = TRUE,
# verbose = TRUE)

# Remove columns with missing values
tmp_rppa_bc_1 <- rppa_bc_1
idx_missing <- which(apply(tmp_rppa_bc_1, 2, function(x) { length(which(is.na(x))) }) > 0)
tmp_rppa_bc_1 <- tmp_rppa_bc_1[, -c(idx_missing)]

data <- tmp_rppa_bc_1
prior <- network_ref 
n_prot <- ncol(proteomic_responses)
k_fold <- 2
proteomic_responses <- proteomic_responses
mab_to_genes <- mab_to_genes
dist_file <- NULL
max_dist <- max_dist
verbose <- verbose
boot_time <- 2
crit_cv <- "BIC"
algorithm <- "data_driven"

a <- sort(colnames(data))
b <- sort(colnames(proteomic_responses))
intersect(a, b)

sort(apply(data, 2, function(x) {
  length(which(is.na(x)))
}))

# Checks: 1) is label in mab? 2) then subset 

x1 <- run_cross_validation(data = data,
                           prior = prior,
                           n_prot = n_prot,
                           k_fold = k_fold,
                           proteomic_responses = proteomic_responses,
                           mab_to_genes = mab_to_genes,
                           dist_file = dist_file,
                           max_dist = max_dist,
                           verbose = verbose,
                           boot_time = boot_time,
                           crit_cv = crit_cv,
                           algorithm = algorithm)

optimize_param <- targetscore::optimize_parameter_dat(
  data = tmp$data_train,
  rho = tmp$rho
)

outersect(colnames(proteomic_responses), colnames(data))
setdiff(colnames(proteomic_responses), mab_to_genes$AntibodyLabel)
setdiff(colnames(data), mab_to_genes$AntibodyLabel)

  
# @param proteomic_responses input drug perturbation data. With columns as antibody, rows as samples.
