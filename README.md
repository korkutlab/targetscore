# Data Formats
zeptosens output file:
sampleno_sample_treatment_dose_time_replicate_date_notes
e.g., S001_skmel133-melanoma_RAFi:MEKi_2000nM:50nM_1000min:1000min_rep1_08272015_this-is-my-note
rules: no space/no punctuation no special characters. Just use letters and numbers. Put all information relating to cell line name, organism name etc. in "sample"" (use "-" as word seperator). Dose in nM (unit “nM"), time in minutes (unit “min”). replicate as "repX” ( rep1 ,rep 2,rep3) .time format month/day/year.

For combos use “:”
time: the drug incubation time.
If delayed experiment time becomes:
t1:t2

_5min:200min_

If you do not have any of these info types, then type “NA”. If you do not have a note
S001_skmel133-melanoma_RAFi:MEKi_2000nM:50nM_1000min:1000min_rep1_08272015_NA
# Installation and Usage 

    install.packages("devtools")
    
    library(devtools)
    install_bitbucket("cbio_mskcc/zeptosensPkg",
        auth_user="discoverUser",
        password="discoverUserPassword",
        build_vignette=TRUE,
        dependencies=TRUE,
        args="--no-multiarch")
        
    library(zeptosensPkg)
