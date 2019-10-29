# Citatation for Data and Algorithm Usage

Please cite the following when using these data: Anil K.et al ...

# Antibody Map File
## File Descriptions

antibody_map_file is the Gene Name list database developed by [Anil Korkut Group](https://odin.mdacc.tmc.edu/~akorkut/#/home), which provided the onsite AnitobodyLabel of local upload files with the Gene_Symbol used within the database.  While antibody_map_file also provided the information of phosphorylation activation/deactivation.

* AntibodyLabel: The AntibodyLabel is the label of local data for each Antibody. The AnitibodyLabel will serve as the label
          alongside all calculation.
* Source: Source of the Antibody. The provided antibody_map_file contains two sources including MDACC as MD Anderson Cancer Center and MSKCC standing for Memorial Sloan Kettering Cancer Center.
* NodeName: MDACC standardized Antibody Name.
* Gene_Symbol: Corresponding HGNC symbol and Ensembl ID
* Sites: phosphorylation site. na stands for no phosphorylation. One site phosphorylation, for example:S473 for Akt_pS473, Two site phosphorylation, for example :Y1234|Y1235 for c.Met_pY1234_Y1235
* Effect: The effect of Phosphorylation. Including:
  * c : no phosphorylation
  * a : activation
  * i : inhibition

## Download

[antibodyMapfile_localfile](data/antibodyMapfile.txt)

# Global signaling data
## File Description

Proteomic datasets, which capture the signaling co-variations serve as  the experimental constraint for network inference. Such datasets can be publicly available (e.g., TCGA data) or custom generated (drug response data, Korkut et al, Elife).The provided example datas

* Columns: HGNC symbol and Ensembl ID
* Rows: Patients Samples

## DownLoad
Here provided an example from Public database TCGA of Breast Cancer.Protein level expression data for all genes, Log2 transformed.

[TCGA_BRCA_localfile](data/TCGA-BRCA-L4_1.csv)

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

[fs_file_localfile](data/Cosmic.txt)
