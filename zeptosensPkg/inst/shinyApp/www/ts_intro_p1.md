# Introduction to TargetScore

## Code Availability 
[Code available as the TargetScore R Package](http://www.git.com)
        
## Why we developed TargetScore?

Targeted therapies have been substantially successful in treatment of diverse cancer types 
_(Gu et al., 2016; Manzano et al., 2016; Mayekar and Bivona, 2017; Ohmoto and Yachida, 2017; Tapia Rico et al., 2017)_.However, resistance to therapy is virtually inevitable and can manifest as a lack of response to therapy (intrinsic)or disease progression after temporary response (acquired resistance) _(Holohan et al., 2013)_.A recurrent mechanism of resistance is activation of compensatory oncogenic pathways (e.g., via feedback loops in short-term or secondary oncogenic alterations in long-term) in response to targeting a genomic aberration _(Holohan et al., 2013; Niederst and Engelman, 2013)_ A relatively simple way to interrogate adaptive responses to targeted agents is to rank changes in mRNA and (phospho)protein expression for individual genes based on high throughput omics data.  However, the rank-based approach cannot capture collective changes that lead to robust phenotypic transitions such as drug resistance.It is also more likely to detect druggable targets within collectively functioning network modules compared to individual genes/proteins. Here, with this motivation, we developed a method to analyze drug response, collective pathway adaptation mechanisms and discover effective drug combinations. **Shown as in Fig.**

## Hypothesis in Target Score Method:

1) To overcome resistance, use combination therapy to target the drug-stimulated pathway.

2) The combination and collective result from pathway may be a better predicter in drug response analysis.

## What does Target Score provided?

Target Score currently have three functional modules which are: Reference network contruction; target score calculation; calculated result visualization. The detailed information of target score calculation process is explained in Help page.