# TargetScore  Tutorial

The tutorial describes running the TargetScore using the package. Users can opt to use the web interface. On the web interface, the data can be uploaded and similar options are provided on dropdown menus.  

## Download and install
For the package and instructions:
https://github.com/korkutlab/targetscore/tree/master/targetscore

# Load the library and input files

Step 1. Load the targetscore library 

Step 2. Definitions and parameters
 Set the maximum pathway neighborhood distance (recommended value: 1)
 Load the antibody Map file to match gene and protein names as well as protein phosphorylation annotations
 
Step 3. Read sample specific drug response data (Log2 normalized to no drug condition). This is needed fo rTargetScopre calculations

Step 4. Read the proteomic data for network inference (e.g., TCG?A RPPA datasets or any proteomic constraints that represent the model system).

Step 5. Prior information from SignedPC or other databases
 
 
Step 6. Functional score file (fs) to annotate proliferative/survival and anti-proflirative/death signals

See sample files for data formats

Tutorial_1.png

# Reference network inference
The reference network model captures the signaling interactions underlying the drug responses.
Use the data from (Step4) and  priors (Step5) to infer the reference network
options
 targetscore::predict_bio_network -> inference based on nly priors from databases (or user set)
 targetscore::predict_dat_network -> inference based on data only with glasso algorithm
 targetscore::predict_hybrid_network -> Inference based on both data and prior constrints (see Wang et al)
 
 **Literature-Based**
Construct the network through Pathway Commons database. Of interest to this study, biochemical reactions, posttranslational modifications, and complex formation, which are all    encoded in BioPAX language are taken into account. [Pathway Commons](https://www.pathwaycommons.org). The pathway interactions that match to the proteomic species profiled by the RPPA data are defined with 4 relation types: phosphorylation, dephosphorylation, expression upregulation, and expression downregulation. The list of gene names and corresponding posttranslational modifications (phosphorylation) is provided to the TargetScore package. The interactions between the molecules of interest are extracted using the signedPC module in BioPAX.The extracted interaction set is evaluated with expert curation and serves as the reference network.


**Graphical LASSO**

Data-driven reference network captures the molecular associations and inter-tumor heterogeneity across a population of samples with shared characteristics (e.g., ovarian cancer patients). Proteomic datasets, which capture the signaling co-variations serve as the experimental constraint for network inference. Such datasets can be publicly available (e.g., TCGA data) or custom generated _(drug response data, Korkut et al, Elife)_. We inferred from glasso algorithm and provided here as the data-driven algorithm in constructing network."),
 [Graphical Lasso algorithm](http://statweb.stanford.edu/~tibs/ftp/glasso-bio.pdf)Glasso generates a partial correlation based network by estimating a sparse inverse of the covariance matrix using an L1 penalty. The algorithm estimates the precision (i.e., inverse covariance) matrix through maximization of the log-likelihood function by applying the L1 penalty"),

 **Hybrid**

We have modified the glasso algorithm and developed the Adjusted-glasso algorithm to infer the reference network model using priors and the directionality information. Adjusted-glasso algorithm introduces the biology prior information from Pathway Commons Database. The introduction of this symmetric prior information matrix can be seen as an adjustment to the inference penalty and allows us to apply varying amounts of penalties to different elements in the precision matrix based on prior interactions. The developed algorithm gives a larger probability for the pre-existing/experimentally validated interactions from Pathway Commons database while taking disease-specific global signaling into account.

 The sif file can be visualised wit network analysis tools (e.g., Cytospcape)
 Tutorial_2.png
 
 # TargetScore calculation 
 (Equation 4 in Wang et al, 2021)
 
Input files and parameters are the proteomic responses, network file (wk), functional scores (fs), nperm sets the number of disstributionss to calculate a null model for statistical assessment.

target score (TS) that quantifies the adaptive pathway responses to a perturbation as a sum of the response from each phosphoprotein level and its pathway neighborhood is calculated for each protein in each sample. The calculation combines the cell type-specific drug response data with the reference network model information. High target score identifies genes involved in adaptive response (e.g., upregulation of RTK mRNA expression by MEK inhibitor via a feedback loop and low target score corresponds to the immediate impact of the drug.") 

**Functional Score for Antibodies**
A functional score of +1 is assigned to proteomic entities representing total level and activating phosphorylations of oncogenes or deactivating phosphorylations of tumor suppressors. Similarly, a functional score of (-1) is assigned to total levels and activating phosphorylations of tumor suppressors and inhibitory phosphorylations of oncoproteins. The shiny app provided a default functional score inferred from Cosmic Database [Cosmic Resources](https://cancer.sanger.ac.uk/cosmic) with a portal for users to upload the self-defined functional score to override.

**TargetScore Calculation method**

Target Score currently provided two methods in the calculation. One listed as Line by Line and the other listed as Pooled.

* **Line by Line**  
Calculation line by line limited the calculation by putting the number of doses as one and calculate target score for each 
every line for the provided Perturbation Response File dataset.

* **Pooled**
Calculation Pooled provided the calculation by putting the number of doses as the number of rows which sum up target score for each 
every line for the provided Perturbation Response File dataset.

**Permutation Number**

To eliminate the connectivity bias, we assess the significance of Target Score for each proteomic entity. For this purpose, the probability of observing a Target Score is calculated over a fixed reference network structure and drug response data with randomized protein labels. The randomized data is generated by random sampling of proteomic responses across all conditions for each antibody. Permutation Number (Default at 25) is the number of randomly bootstrapped data sets, and the null distribution of Target Scores for a given network topology is calculated.

**Maximum Network Distance**

Maximum Network Distance between two nodes. Which limits the distance between nodes from the reference network that will be included in the Target Score calculation.

## TargetScore Visualization

We provided three modules for visualization including Network Edgelist from reference network construction, Heatmap for calculated target score, and Target Score Volcano Plot.

 
 Tutorial_3.png
 
 # Output files
 TS_sample.csv -> The target scores for each protein profiled
 TS_sample_pvalue.csv -> The p-values for each TargetScaore against a null model generated by bootstrapping with label randomized data
 TS_sample_qvalue.csv -> The FDR-adjussted (BH-method) q-values based on the p-values

Tutorial_4.png




