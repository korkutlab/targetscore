# Citation for Data and Algorithm Usage

Please cite the following when using the example data/source code/TargetScore calculator: 
Anil Korkut, Heping Wang, Augustin Luna.

# Target Score Result File

## File Descriptions

Target Score Result File is the target score calculated result provided by the user. The data file should be an RDS data file as a list with ts_r: the target score calculated result list, proteomic_responses: Perturbation Response, fs_dat:Functional Score, mab_to_genes: Mapping Gene list(Example:Mapping File), network: network used by the calculation. The Input Panel is provided for users who are insterested in Shiny App visualization part and want to short the calculation time.


# Perturbation Response File

## File Descriptions

Perturbation Response File is the drug perturbation data provided by the user. We have provided an example drug perturbation response file here which provided by Anil Korkut lab. The data file contains Antibody Label as column and drug perturbation treatments as rows.

## Download
Here provided an example from Anil Korkut lab of cell line BT474 (Breast Cancer cell line).Protein level expression data for ~300 proteins, Log2 transformed and median centered.
[Perturbation Response file directory](https://github.com/HepingWang/zeptosensPkg/blob/master/zeptosensPkg/inst/test_data/BT474.csv)

# Mapping File

## File Descriptions

Mapping file is the Gene Name list database developed by [Anil Korkut Group](https://odin.mdacc.tmc.edu/~akorkut/#/home), which provided the onsite AnitobodyLabel of local upload files with the Gene_Symbol used within the database.  While antibody_map_file also provided the information of phosphorylation activation/deactivation.

* AntibodyLabel: The AntibodyLabel is the label of local data for each Antibody. The AnitibodyLabel will serve as the label alongside all calculation.
* Source: Source of the Antibody. The provided antibody_map_file contains two sources including MDACC as MD Anderson Cancer Center and MSKCC standing for Memorial Sloan Kettering Cancer Center.
* NodeName: MDACC standardized Antibody Name.
* Gene_Symbol: Corresponding HGNC symbol and Ensembl ID
* Sites: phosphorylation site. na stands for no phosphorylation. One site phosphorylation, for example:S473 for Akt_pS473, Two site phosphorylation, for example :Y1234|Y1235 for c.Met_pY1234_Y1235
* Effect: The effect of Phosphorylation. Including:
  * c : no phosphorylation
  * a : activation
  * i : inhibition

## Download
Here provided an example file of Mapping file which used by Korkut lab currently. 
[Mapping file directory](https://github.com/HepingWang/zeptosensPkg/blob/master/zeptosensPkg/inst/targetScoreData/antibodyMapfile_08272020.csv)

# Background Network File

## File Description
Background Network File, which capture the signaling co-variations serve as  the experimental constraint for network inference. Such datasets can be publicly available (e.g., TCGA data) or custom generated (drug response data, Korkut et al, Elife).The provided example datas

* Columns: HGNC symbol and Ensembl ID
* Rows: Patients Samples

## DownLoad
Here provided an example from Public database TCGA of Breast Cancer.Protein level expression data for all genes, Log2 transformed.

[TCGA BRCA file directory](https://github.com/HepingWang/zeptosensPkg/blob/master/zeptosensPkg/inst/test_data/TCGA-BRCA-L4.csv)

# Functional Score file
## File Description

A functional role is assigned as a numeric score to proteomic entities.

* there is evidence for an entity’s function as an oncogene or tumor suppressor in cancer.
* A central basis for this cancer role comes the curated and widely-used COSMIC database.
* Two different functional score value were assigned: oncogene as +1,and tumor suppressor as -1.
* it is also possible for users to manually alter these scores referring from literature or through expert editing, if necessary."),
* gene : HGNC symbol and Ensembl ID
* fs : Corresponding functional score

## Download
Here provided an example file of functional score file from the Catalogue Of Somatic Mutations In Cancer project resource where corresponding functional score of listed proteins:
0 for dual or unknown proteins;
-1 for tumor suppressors;
1 for onco-proteins.
[Functional Score file directory](https://github.com/HepingWang/zeptosensPkg/blob/master/zeptosensPkg/inst/extdata/Cosmic.txt)
