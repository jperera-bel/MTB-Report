# MTB-Report

The user can generate Molecular Tumor Board (MTB) reports for TCGA samples. MTB reports include a filtered list of actionable variants. To do so, it follows the following steps:

1. Input SNVs, CNVs and fusions are queried against databases of actionable variants:
    - Dienstmann et al., Cancer discovery 5.2 (2015), v19
    - Griffith et al., Nat Genet (2017), version 1 July 2017
    - Van Allen et al., Nat Med (2014), v3
    - Meric-Bernstam et al., J Natl Cancer Inst. (2015)
2. Matching variants are then classified into levels of evidence
3. Finally, a pdf report is generated

## Dependencies
```r
install.packages("knitr")
install.packages("stringr")
install.packages("xtable")
install.packages("ggplot2")
install.packages("pander")
install.packages("timeSeries")
devtools::install_github("mariodeng/FirebrowseR")
```

Requires LaTeX

## Files Description

#### script.r 
Main script. Allows the user to filter SNVs and CNVs using gene-drug public 
databases. Then classifies the variants into levels of evidence 
and finally generates report in pdf format. 


#### /data folder
Contains the databases used to filter for actionable variants
and an in-house file with cancer type equivalences between the databases


#### /helpers folder
Contains three files with the helper functions used by script.r

