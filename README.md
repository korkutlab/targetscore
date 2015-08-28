# Data Formats
# Format 

sampleNumber_sample_treatment_dose_time_replicate_date_notes

# Complete Example
S001_skmel133-melanoma_RAFi:MEKi_2000nM:50nM_1000min:1000min_rep1_20150827_this-is-my-note

# Sample Rules

Put all information relating to cell line name, organism name etc in "sample" 

* No spaces
* No punctuation 
* No special characters 

* Use only use numbers [0-9], letters [A-Za-z]
* Use dashes [-] between words

## Unit Rule 

* Dose in nM (unit "nM")
* Time in minutes (unit "min"); time is the drug incubation time 
** Example: 10min, 360min (6 hours), 1440min (24 hours) 

## Date Format Rule 

* Year, Month, Day: YYYYMMDD
** Example: 20150828

## Replicate Rule 

* Replicate as "repX" 
** Example: rep1 ,rep2, rep3

## Drug Combination Rule

* Use colon [:] between multiple drugs, times, doses 
* Provide information for each drug even if the same
** Example: RAFi:MEKi, 2000nM:50nM, 1000min:1000min

### Delayed Experiment Time 

* Delayed experiment time becomes: t1:t2
** Example: 5min:200min

## Missing Information

* If you do not have any of these info types, then type "NA". 
* Example without note: S001_skmel133-melanoma_RAFi:MEKi_2000nM:50nM_1000min:1000min_rep1_20150827_NA

# Total Protein Levels Format 
Sample Name
Protein Concentration 

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
