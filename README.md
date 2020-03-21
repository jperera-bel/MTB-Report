# MTB-Report

The user can generate Molecular Tumor Board (MTB) reports for TCGA samples. MTB reports include a filtered list of actionable variants. To do so, it follows the following steps:

1. Input SNVs, CNVs and fusions are queried against databases of actionable variants:
    - [Gene Drug Knowledge Database](https://www.synapse.org/#!Synapse:syn2370773): Dienstmann et al., Cancer discovery 5.2 (2015), v20.0
    - [CIViC](https://civic.genome.wustl.edu/): Griffith et al., Nat Genet (2017), release 01-Mar-2020
    - [TARGET](http://archive.broadinstitute.org/cancer/cga/target): Van Allen et al., Nat Med (2014), v3
    - Meric-Bernstam et al., J Natl Cancer Inst. (2015)
2. Matching variants are then classified into levels of evidence
3. Finally, a pdf report is generated

If you use MTB report, please cite this pulication:

<b>Perera-Bel J</b>, Hutter B, Heining C, Bleckmann A, Fröhlich M, Fröhling S, Glimm H, Brors B, Beißbarth T. <i>From somatic variants towards precision oncology: Evidence-driven reporting of treatment options in molecular tumor boards</i>. Genome Med. 2008 15;10(1):18. <a target="_blank" href="https://doi.org/10.1186/s13073-018-0529-2">doi: 10.1186/s13073-018-0529-2.</a>

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

Requires LaTeX an Texinfo. In Linux, install them with:
```
sudo apt-get install texlive-full
sudo apt-get install texinfo
```
## Usage

Open the main R script (script.r) in R or RStudio. The user can change the default patient ID (TCGA ID) to any TCGA sample which is available through Firebrowse. 

This script can also be used to analyze user-defined data as long as it follows the same data structure:

#### SNVs

| Hugo_Symbol   | Variant_Classification| Protein_Change  |
| ------------- |:----------------------:| --------------:|
| EGFR          | missense mutation      | T790M          |
| APC           | nonsense mutation      |   K670*        |

#### CNVs

| Hugo_Symbol   | CN alteration|
| ------------- |:------------:|
| HER2          | amplification| 
| FBN3          | deletion     | 

*Gene names* must be Hugo Symbols. 

*Variant Classification* comprises the following levels: {Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region}. 

*Protein Chanage* must be e.g. T790M

*CN alteration* must be one of the foollowing levels: {amplification, deletion}_

## Files Description

#### script.r 
Main script. Allows the user to filter SNVs and CNVs using gene-drug public 
databases. Then classifies the variants into levels of evidence 
and finally generates report in pdf format. 


#### /data folder
Contains the databases used to filter for actionable variants
and an in-house file with cancer type equivalences between the databases


#### /helpers folder
Contains four files with the helper functions used by script.r

