# TargetScore Tutorial

The tutorial describes running the TargetScore using the package. Users can opt to use the web interface. On the web interface, the data can be uploaded and similar options are provided on dropdown menus.  

# Download and Install

For the package and instructions:

https://github.com/korkutlab/targetscore/tree/master/targetscore

# Steps to Load the Library and Input Files

1. Load the targetscore library 
2. Definitions and parameters
  * Set the maximum pathway neighborhood distance (recommended value: 1)
  * Load the antibody Map file to match gene and protein names as well as protein phosphorylation annotations
  * Antibody Map File: FIXME link to the antibody map file
  * File Description: FIXME link to file descriptions
3. Read sample specific drug response data (log2-normalized to no drug condition). This is needed fo targetScore calculations
  * Sample File: link to sample file
  * File Description: link to file descriptions
4. Read the proteomic data for network inference (e.g., TCGA RPPA datasets or any proteomic constraints that represent the model system).
  * Sample File: link to sample file
  * File Descriptions: link to file descriptions
5. Prior information from SignedPC or other databases
  * Signed PC File: link to signed PC
  * File Description: link to file descriptions
6. Functional score file (fs) to annotate proliferative/survival and anti-proflirative/death signals
  * Functional Score File: link to fs.txt file
  * File Description: link to file descriptions

## Code
```
#LOAD THE PACKAGE/LOAD DATA AND INPUT PARAMETERS
#1. Load Library
library(zeptosensPkg)
# 2. Set the maximum pathway neighborhood distance in TS calculations)
maxDist <- 1 
# Read in Protein list file
mab_to_genes <- read.table(system.file("targetscoreData", "antibodyMapfile.csv", package = "zeptosensPkg"),
                           sep = ",",
                           header = TRUE,
                           stringsAsFactors = FALSE)

#3. Read proteomic response for cellline1 (This is needed in the nertwork inference to import protein names & as data for TS calculation)
x_1 <- read.csv( "data_TS/HCC1954.csv",row.names = 1)  #nProt=304


#4. Read The Global Signaling file for BRCA (Source: TCGA, data for reference network inference)
RefNetworkData <- read.csv( "data_TS/TCGA-BRCA-L4.csv",row.names = 1)

# 5. Read priors (Source: SignedPC, Babur et al, 2019)
prior <- read.csv(file="data_TS/Bio_network.csv",row.names = 1)
##

# 6. Get the functional scores
fs_override_org <- readRDS(system.file("test_data_files", "fs_value_file.rds",
                                       package = "zeptosensPkg"
))

fs_value <- zeptosensPkg::get_fs_vals(
  n_prot = ncol(x_1), proteomic_responses = x_1,
  mab_to_genes = mab_to_genes, fs_override = fs_override_org
)
```

# Reference Network Inference
The reference network model captures the signaling interactions underlying the drug responses.
Use the data from (Step 4) and  priors (Step 5) to infer the reference network
options
* targetscore::predict_bio_network -> inference based on only priors from databases (or user set)
* targetscore::predict_dat_network -> inference based on data only with glasso algorithm
* targetscore::predict_hybrid_network -> Inference based on both data and prior constrints (see Wang/Luna et al)
 
## Literature-Based
Construct the network through Pathway Commons database. Of interest to this study, biochemical reactions, posttranslational modifications, and complex formation, which are all    encoded in BioPAX language are taken into account. [Pathway Commons](https://www.pathwaycommons.org). The pathway interactions that match to the proteomic species profiled by the RPPA data are defined with 4 relation types: phosphorylation, dephosphorylation, expression upregulation, and expression downregulation. The list of gene names and corresponding posttranslational modifications (phosphorylation) is provided to the TargetScore package. The interactions between the molecules of interest are extracted using the signedPC module in BioPAX.The extracted interaction set is evaluated with expert curation and serves as the reference network.

## Graphical LASSO
Data-driven reference network captures the molecular associations and inter-tumor heterogeneity across a population of samples with shared characteristics (e.g., ovarian cancer patients). Proteomic datasets, which capture the signaling co-variations serve as the experimental constraint for network inference. Such datasets can be publicly available (e.g., TCGA data) or custom generated _(drug response data, Korkut et al, Elife)_. We inferred from glasso algorithm and provided here as the data-driven algorithm in constructing network."), [Graphical Lasso algorithm](http://statweb.stanford.edu/~tibs/ftp/glasso-bio.pdf)Glasso generates a partial correlation based network by estimating a sparse inverse of the covariance matrix using an L1 penalty. The algorithm estimates the precision (i.e., inverse covariance) matrix through maximization of the log-likelihood function by applying the L1 penalty"),

## Hybrid
We have modified the glasso algorithm and developed the prior-glasso algorithm to infer the reference network model using priors and the directionality information. The prior-glasso algorithm introduces the biologically relavant prior information from Pathway Commons Database.  

### Code
```
#REFERENCE NETWORK INFERENCE
#prior-glasso inference
network.Hyb<-zeptosensPkg::predict_hybrid_network(
  n_prot = 304 ,
 proteomic_responses = x_1,
  prior = prior,
  data =tcga_data,
  cut_off = 0.05
)
#write the network matrix
write.csv(network.Hyb$wk,file="RefNetwork.csv")

hybnet<-read.csv(file="RefNetwork.csv",row.names = 1)
#write the network SIF file
hybnet_sif<-zeptosensPkg::create_sif_from_matrix(t_net = hybnet,col_genelist = colnames(hybnet),row_genelist = rownames(hybnet))
write.table(hybnet_sif,file="ReferenceNetwork_sif.txt",quote = F,sep = "\t",row.names = F)
```
 
# TargetScore Calculation 
Input files and parameters are the proteomic responses, network file (wk), functional scores (fs), nperm sets the number of disstributionss to calculate a null model for statistical assessment.

TargetScore (TS)  quantifies the adaptive pathway responses to a perturbation as a sum of the response from each phosphoprotein level and its pathway neighborhood is calculated for each protein in each sample. The calculation combines the cell type-specific drug response data with the reference network model information. High TargetScore identifies processes involved in adaptive response (e.g., upregulation of RTK mRNA expression by MEK inhibitor via a feedback loop and low target score corresponds to the immediate impact of the drug.) 

## Functional Score for proteins
A functional score of +1 is assigned to proteomic entities representing total level and activating phosphorylations of oncogenes or deactivating phosphorylations of tumor suppressors. Similarly, a functional score of (-1) is assigned to total levels and activating phosphorylations of tumor suppressors and inhibitory phosphorylations of oncoproteins. The shiny app provided a default functional score inferred from Cosmic Database [Cosmic Resources](https://cancer.sanger.ac.uk/cosmic) with a portal for users to upload the self-defined functional score to override.
 
## Permutation Number
To eliminate the connectivity bias, we assess the significance of Target Score for each proteomic entity. For this purpose, the probability of observing a Target Score is calculated over a fixed reference network structure and drug response data with randomized protein labels. The randomized data is generated by random sampling of proteomic responses across all conditions for each antibody. Permutation Number (Default at 25) is the number of randomly bootstrapped data sets, and the null distribution of Target Scores for a given network topology is calculated.

### Code
```
#TARGET SCORE CALCULATION
proteomic_responses_1 <- x_1
length(proteomic_responses_1)

# Calculate Target Score
ts <- array(0, dim = c(dim(proteomic_responses_1)[1], dim(proteomic_responses_1)[2]))
ts_p <- array(0, dim = c(dim(proteomic_responses_1)[1], dim(proteomic_responses_1)[2]))
ts_q <- array(0, dim = c(dim(proteomic_responses_1)[1], dim(proteomic_responses_1)[2]))

for (i in seq_len(dim(proteomic_responses_1)[1])) {
  results <- zeptosensPkg::get_target_score(
    wk = network.Hyb$wk,
    wks = network.Hyb$wks,
    dist_ind = network.Hyb$dist_ind,
    inter = network.Hyb$inter,
    n_dose = 1,
    n_prot = dim(proteomic_responses_1)[2],
    proteomic_responses = proteomic_responses_1[i, ],
    n_perm = 1000,
    verbose = FALSE,
    fs_dat = fs_value
  )
  ts[i, ] <- results$ts
  ts_p[i, ] <- results$pts
  ts_q[i, ] <- results$q
}
```
 
# Output Files
* TS_sample.csv: The target scores for each protein profiled
* TS_sample_pvalue.csv: The p-values for each TargetScore against a null model generated by bootstrapping with label randomized data
* TS_sample_qvalue.csv: The FDR-adjusted (BH-method) q-values based on the p-values

## Code
```
#OUTPUT FILES
colnames(ts) <- colnames(proteomic_responses_1)
ts <- data.frame(rownames(proteomic_responses_1), ts)
write.csv(ts,file="data_TS/TS_HCC1954.csv",row.names = F)

colnames(ts_p) <- colnames(proteomic_responses_1)
ts_p <- data.frame(rownames(proteomic_responses_1), ts_p)
write.csv(ts_p,file="data_TS/TS_HCC1954_pvalue.csv",row.names = F)


colnames(ts_q) <- colnames(proteomic_responses_1)
ts_q <- data.frame(rownames(proteomic_responses_1), ts_q)
write.csv(ts_q,file="data_TS/TS_HCC1954_qvalue.csv",row.names = F)
```

# TargetScore Visualization
We provided three modules for visualization including network edgelist from reference network construction, heatmap for calculated target score, and Target Score volcano plot.
