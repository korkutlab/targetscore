# TargetScore 
TargeScore is a statiscal network modeling and analysis algorithm

The algorithm 
(i) reveals and visualizes collective/pathway activity changes involved in adaptive resistance as drug-activated network modules of signaling;  
(ii) nominates combination therapies to down-regulate the resistance pathways. 

This the TargetScore web application and enables calculation of drug-activated network modules of signaling. 

The minimal molecular data requirement is molecular responses to a single sample treated with a single agent. 

The comparison of TargetScore values over different cancer samples (e.g., sensitive vs. resistant), across different drug doses or time points reveals how drugs rewire pathway activities in time, dose and sample space.

TS_figure.png

The algorithm involves the following steps: 1. Molecular profiling of the cellular response to a perturbation; 2. Inference of a reference network (e.g., breast cancer signaling network) that captures potential relations between measured proteomic entities across diverse samples in a specific class; 3. Quantification of a sample and context-specific adaptation score (i.e., a TargetScore value) which links protein interactions to drug response on the reference network using molecular drug response data. The TargetScore values reflect adaptive responses and are calculated for each protein under each condition. It is quantified as the network interaction weighted sum of the "self-change" of the corresponding entity and the change in the pathway/network neighborhood in response to targeted perturbations. 4. Identification of network modules (collective changes) that have a significantly high TargetScore (i.e., collectively participate in adaptive responses) in each sample. The network modules involved in adaptive responses are determined by mapping the TargetScore values back on the reference network and extracting the connected sub-networks enriched with high TargetScore values. 5. Selection of actionable targets that participate in adaptive responses in a given sample and test drug combinations in pre-clinical models. 

#Code Availability
Code available as the TargetScore R Package (github link)

# Feedback

We appreciate any feedback/suggestions you may have. Please forward feedback to [Anil Korkut](mailto:akorkut@mdanderson.org)
