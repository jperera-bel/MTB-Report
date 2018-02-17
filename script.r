#==============================================================#
#       GENERATION OF MOLECULAR TUMOR BOARD (MTB) REPORT       #
#--------------------------------------------------------------#
#   This script filters SNVs and CNVs using gene-drug public   #
#   databases. Then classifies the variants into levels of evi-#
#   dence and finally presents the results in a pdf report     #
#==============================================================#


## Dependencies
library(knitr)
library(stringr)
library(xtable)
library(pander)
library(ggplot2)
library(timeSeries)
library(FirebrowseR)

######################################################
## Load everything necessary for downstream analysis #
######################################################

## Load helper functions
source("helpers/get_druggable.r")
source("helpers/get_levels.r")

# Cancer type synonyms between databases
synonyms = read.csv('data/cancer_types.csv', header = TRUE,sep="\t")

###################################################
## Get clinical, SNVs and CNVs from a TCGA sample #
###################################################

# As an example, in this script we will generate a report for a TCGA sample
# The user can generate reports for other TCGA samples by changing the ID
# User defined data can be analyzed as long as it follows the same data structure

patient  = "TCGA-34-5239" # lusc

# Here some other example samples to try
#patient  = "TCGA-AA-A03F" # coad with KRAS mut
#patient  = "TCGA-AA-A01S" # coad with KRAS wt
#patient  = "TCGA-B0-4827" # kirc
#patient  = "TCGA-CQ-5325" # hnsc
#patient  = "TCGA-AR-A1AR" # brca
#patient  = "TCGA-D1-A0ZZ" # ucec


clinical = Samples.Clinical(format = "csv", tcga_participant_barcode = patient)
SNV      = setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("Hugo_Symbol","Variant_Classification","Protein_Change"))
CNV      = setNames(data.frame(matrix(ncol = 2, nrow = 0)),c("gene","cn_alteration"))
SNV      = Analyses.Mutation.MAF(format = "csv",tcga_participant_barcode = patient)
CNV      = Analyses.CopyNumber.Genes.Thresholded(format = "csv",tcga_participant_barcode = patient,page_size =25000) # all genes

if( nrow(SNV) == 0 & nrow(CNV) == 0) warning('The patient has no SNVs nor CNVs. Not worth continuing with the analysis!')

## Get cancer type synonyms in the databases
cancer           = clinical$cohort # Some samples belong to two cohors! Check!
cancer_GDKD      = unique(as.character(synonyms[grep(cancer,synonyms$tcga_cancer,ignore.case = T),"knowledge"]))
cancer_CIVIC     = sapply(cancer,function(x) as.character(na.omit(synonyms[grep(x,synonyms$tcga_cancer,ignore.case = T),"civic"])[1]))
cancer_CIVIC     = paste(cancer_CIVIC,collapse = ",")
cancer_CIVIC     = unique(strsplit(cancer_CIVIC,",")[[1]])

#####################################
## Filter SNVs and CNVs by database #
#####################################

## Prepare inputs
# Important that SNVs variable has 3 columns (Gene name, Variant type, aa change)
SNV                = SNV[SNV$Variant_Classification!='Silent',]
SNV$Protein_Change = gsub("p.","",SNV[,"Protein_Change"])
SNV                = SNV[,c("Hugo_Symbol","Variant_Classification","Protein_Change")]
SNV                = SNV[which(SNV[,"Protein_Change"] != ""),]
SNV                = unique(SNV)

# Important that CNVs variable has 2 columns (Gene name, CNV type)
CNV = CNV[,c("gene","cn_alteration")]
CNV[CNV[,2] ==-1,2] = 0 # Here we don't consider broad/low level deletions
CNV[CNV[,2] == 1,2] = 0 # Here we don't consider broad/low level aplifications
CNV[CNV[,2] ==-2,2] = "deletion"
CNV[CNV[,2] == 2,2] = "amplification"
CNV                 = CNV[which(CNV[,2] != 0),]
CNV                 = unique(CNV)

if( nrow(SNV) == 0 & nrow(CNV) == 0 ) warning('The patient has no SNVs nor CNVs. Not worth continuing with the analysis!')

#### GDKD DB
druggableGDKD = data.frame()
druggableGDKD = match_SNV_GDKD(SNV)
druggableGDKD = rbind(druggableGDKD,match_CNV_GDKD(CNV))
druggableGDKD = rbind(druggableGDKD,match_WT_GDKD(SNV,CNV,cancer_GDKD))
rownames(druggableGDKD) = NULL

#### CIVIC
druggableCIVIC = data.frame()
druggableCIVIC           = match_SNV_CIVIC(SNV)
druggableCIVIC           = unique(rbind(druggableCIVIC,match_CNV_CIVIC(CNV)))
rownames(druggableCIVIC) = NULL


#### TARGET DB
druggableTARGET = data.frame()
druggableTARGET = match_TARGET_MERIC(SNV,CNV)

if( nrow(druggableGDKD) == 0 & nrow(druggableCIVIC) == 0 & nrow(druggableTARGET)==0) warning('No gene-drug interations were found. Not worth continuing with the analysis!')


#########################################################
## Classify filtered variants by Levels of Evidence    ##
#########################################################

### KNWOLEDGE
levelsGDKD = c()
levelsGDKD = get_levels_GDKD(druggableGDKD,cancer_GDKD)
#print(levelsGDKD)

### CIVIC
levelsCIVIC = c()
levelsCIVIC = get_levels_CIVIC(druggableCIVIC,cancer_CIVIC)
#print(levelsCIVIC)

levels = merge_levels(levelsGDKD,levelsCIVIC)

# Homogeneize/clean the final table
table = clean_levels(levels,synonyms,sort_by="drug_freq")

#################################################
## Generate report from table and patient data ##
#################################################

# Rename variables
# Important to use the variable names below
patient = patient
pat     = clinical[1,]
cancer  = cancer[1]
maf     = SNV
cnv     = CNV
A       = table
druggableTARGET = druggableTARGET


# Create report
Sweave2knitr("helpers/Report_tcga.Rnw")
knit2pdf("helpers/Report_tcga-knitr.Rnw", clean = TRUE)

# Now, the pdf report should be created in the same directory as script.R




