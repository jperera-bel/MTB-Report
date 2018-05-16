library(timeSeries)
library(xlsx)

## GDKD: Knowledge_database_v20.0.xlsx" dowloaded from https://www.synapse.org/#!Synapse:syn2370773/files/
## CIViC: 01-May-2018-ClinicalEvidenceSummaries.tsv" dowloaded from https://civic.genome.wustl.edu/releases



	
gdkd         = read.xlsx2('Knowledge_database_v20.0.xlsx', header=TRUE, sheetIndex = 1,colIndex = c(1:45))
civic        = read.delim('01-May-2018-ClinicalEvidenceSummaries.tsv', header=T,stringsAsFactors = F,sep="\t",quote = "")
#target_meric = read.xlsx2('TARGET_db_v3_02142015.xlsx', header=TRUE, sheetIndex = 1)

# Small changes
gdkd$Gene  = gsub(" ", "", gdkd$Gene)
for (c in colnames(gdkd)){
	gdkd[,c]  = gsub("^ | $", "", gdkd[,c])
}


for (c in grep("PMID",colnames(gdkd))){
  gdkd[,c] = gsub("\\^|_","",gdkd[,c])
}
civic <- civic[which(civic$evidence_type=="Predictive"),]

# Aggregate databases (join variants with same evidence of: disease, association, drug, evidence level)
gdkd_agg   = aggregate(gdkd, by = list(gdkd$Disease,gdkd$Gene,gdkd$Description,gdkd$Association_1,gdkd$Therapeutic.context_1),
                           FUN = function(X) paste(unique(X), collapse=", "))[,6:50]
civic_agg  = aggregate(civic, by = list(civic$gene,civic$disease,civic$drugs,civic$evidence_level,civic$clinical_significance),
                        FUN = function(X) paste(unique(X), collapse=", "))[,6:41]

gdkd  = gdkd_agg
civic = civic_agg


write.table(gdkd, file="GDKD.csv",col.names=T,row.names=F,sep="\t")
write.table(civic, file="CIViC.csv",col.names=T,row.names=F,sep="\t")


