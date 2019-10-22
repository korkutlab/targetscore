## What does TargetScore do?

The algorithm determines the network level, collective responses to perturbations  with a multi-step approach.Target Score Basically have two steps.

**Step1. Contructing Reference Network**
In the first step, we generate a reference network model that capture the signaling interactions underlying the drug responses. The reference network is inferred using either a knowledge based approach with pathway information from multiple signaling databases or a data-driven approach that uses Hybridh Adjusted-Glasso Algorithm (HAGA).There are three provided ways of Contructing Reference Network.

**predictBionetwork Module:**

Contruct the network through Pathway Commons database.Of interest to this study, biochemical reactions, posttranslational modifications, and complex formation, which are all encoded in BioPAX language are taken into account. [Pathway Commons](https://www.pathwaycommons.org).The pathway interactions that match to the proteomic species profiled by the RPPA data are defined with 4 relation types: phosphorylation, dephosphorylation, expression upregulation, and expression downregulation. The list of gene names and corresponding posttranslational modifications (phosphorylation) is provided to the TargetScore package. The interactions between the molecules of interest are extracted using the signedPC module in BioPAX.The extracted interaction set is evaluated with expert curation and serve as the reference network.

**predictDatnetwork Module**

Data driven reference network captures the molecular associations and inter-tumor heterogeneity across a population of samples with shared characteristics (e.g., ovarian cancer patients).Proteomic datasets, which capture the signaling co-variations serve as  the experimental constraint for network inference.Such datasets can be publicly available (e.g., TCGA data) or custom generated _(drug response data, Korkut et al, Elife)_.We inferred from glasso algorithm and provided here as the data driven algorithm in constructing network."),
 [Graphical Lasso algorithm](http://statweb.stanford.edu/~tibs/ftp/glasso-bio.pdf)Glasso generates a partial correlation based network by estimating a sparse inverse of the covariance matrix using an L1 penalty.The algorithm estimates the precision (i.e., inverse covariance) matrix through maximization of the log-likelihood function by applying the L1 penalty"),

 **predictHybnetwork**

We have modified the glasso algorithm and developed the Adjusted-glasso algorithm to infer the reference network model using priors and with directionality information. Adjusted-glasso algorithm introduces the biology prior information from Pathway Commons Database.The introduction of this symmetric prior information matrix can be seen as an adjustment to the inference penalty and allows us to apply varying amounts of penalties to different elements in precision matrix based on prior interactions.The developed algoritm give larger probability for the pre-existing/experimentally validated interactions from Pathway Commons database while take disease specific global signaling into account.

**Step2. Calculating Target Score Through Reference Network**

target score (TS) that quantifies the adaptive pathway responses to a perturbation as a sum of the response from each individual phosphoprotein level and its pathway neighborhood is calculated for each protein in each sample.The calculation combines the cell type-specific drug response data with the reference network model information. High target score identifies genes involved in adaptive response (e.g., upregulation of RTK mRNA expression by MEK inhibitor via a feedback loop _[[CITE]])_ and low target score corresponds to the immediate impact of the drug.") # tags$a(href="http://www.git.com","TargetScore Pakage Bio"),