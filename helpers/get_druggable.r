##============================================================================##
##                     MATCHING GENE-DRUG INTERACTIONS                        ##
##----------------------------------------------------------------------------##
##                                                                            ##
## Set of functions that search for matching variants (SNVs, CNVs and fusion  ##
## genes) in specialized  gene-drug interaction databases, specifically:      ##
##                                                                            ##
## 1) GDKD: Dienstmann et al., Cancer discovery 5.2 (2015): 118-123, v19      ##
## 2) CIVic: Griffith et al., bioRxiv (2016): 072892, downloaded 01/06/2017   ##
## 3) TARGET: Van Allen et al., Nature medicine 20.6 (2014): 682-688, v3      ##
## 4) Meric-Bernstam et al., J Natl Cancer Inst. 107(7) (2015)                ##
##                                                                            ##
##============================================================================##

#############
## 1) GDKD ##
#############


match_SNV_GDKD = function(snv,db = read.delim("data/GDKD.csv",sep="\t")){
  ## Given a list of SNVs, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
  ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## output: returns those Gene Drug Knowledge Database rows matching to the input SNVs.
  
  gdkd=db
  #gdkd = read.delim("data/GDKD.csv",sep="\t")
  gdkd$Patient_variant = ""
  druggable                 = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
  colnames(druggable)       = colnames(gdkd)
  colnames(druggable)       = colnames(gdkd)
  
  ## We apply a narrowing proceedure
  ## Starting with all variants of the database
  ## We filter out genes in the database not mutated in the patient
  gdkd$Gene= gsub(" ", "",gdkd$Gene)
  index               = as.character(gdkd$Gene) %in% as.character(snv[,1])  #Database genes mutated in patient
  druggable           = subset(gdkd,index)                                  #Subset from complete DB
  rownames(druggable) = NULL

  ## From the subset previously filtered by gene
  ## We annotate only matching variants or repurposed ones (NEW)
  if (nrow(druggable) != 0){
  
    for (i in 1:nrow(druggable)){
    
      gene       =  druggable[i,"Gene"]
      p_var_type = as.character(snv[which(as.character(snv[,1])==gene),2]) # Patient's variant type on current gene
      p_variants = unique(as.character(snv[which(as.character(snv[,1])==gene),3])) # Patient's protein change on the current gene
      for (p_variant in p_variants){
      
        ### Match Missense mutations with variant not specified in the database
        if (grepl("missense|mutation|any",druggable[i,"Description"], ignore.case = TRUE) && 
          grepl("mut|any|unknown",druggable[i,"Variant"], ignore.case = TRUE) && 
          grepl("missense|frame_shift",p_var_type,ignore.case = TRUE) ){
          druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
        }
        
        ### Match Missense mutations with specific variant specified in the database
        if (grepl("missense",druggable[i,"Description"], ignore.case = TRUE) && grepl("missense",p_var_type,ignore.case = TRUE) &&
        grepl("mut|any|unknown",druggable[i,"Variant"],ignore.case = TRUE, )==FALSE){
          if (grepl(p_variant,druggable[i,"Variant"])){ 
            druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant) }

          ## Repurposing - New variants         
          else { 
            if (grepl(gene,"BRAF,KIT,FGFR3,KRAS,PDGFRA,EGFR,AKT1,MTOR,ALK",ignore.case = TRUE)){ 
              if (gene=="EGFR" &&  druggable[i,"Variant"] != "T790M" && druggable[i,"Variant"] != "C797S"){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
              if (gene=="MTOR" && druggable[i,"Variant"] != "F2108L"){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
              if (gene=="AKT1" && druggable[i,"Variant"] != "Q79K"){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
              if (gene=="BRAF" && grepl("V600",druggable[i,"Variant"])==F){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
            }
            else {
              druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ") 
            }
          }
        }
        
        ### Match indels, small deletions and insertions - repurposing rule (position is not checked)
        if (grepl("indel|deletion mutation|insertion mutation",druggable[i,"Description"], ignore.case = TRUE) && 
          grepl("ins|del",p_var_type,ignore.case = TRUE) ){
          druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
        }
        
        ### Match Nonsense mut - repurposing rule (deletion, LoF)
        if (grepl("loss-of-function",druggable[i,"Effect"], ignore.case = TRUE) && 
          grepl("nonsense",p_var_type,ignore.case = TRUE) ){
          druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
        }
        
        ### Find matches in variants annotated at exon level (exon X p.XX-XX)
        if (grepl("-",druggable[i,"Variant"])){
          aa_pos = as.numeric(str_match(p_variant,"[0-9]+"))
          range  = as.character(str_match_all(as.character(gdkd[i,"Variant"]),"[0-9]+-[0-9]+")[[1]])
          range  = strsplit(range,"-")
          if (length(range) != 0){
            for (x in 1:length(range)){
              if (aa_pos >= range[[x]][1] && aa_pos <= range[[x]][2]){druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant) } 
            }
          }
        }
      }
    }
  }

  ## Keep all those entries with "Patient_variant"
  index= which(druggable$Patient_variant != "")
  druggable=druggable[index,]   #Subset from first selection
                                       
  return(druggable)
  
}



match_CNV_GDKD = function(cnv,db = read.delim("data/GDKD.csv",sep="\t")){
  ## Given a list of CNVs, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## output: returns those Gene Drug Knowledge Database rows matching to the input CNVs.
  
  gdkd = db
  gdkd$Patient_variant = ""
  druggable                 = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
  colnames(druggable)       = colnames(gdkd)
  
  
  ## We apply a narrowing proceedure
  ## Starting with all variants of the database
  ## We filter out genes in the database with no CNVs in the patient
  index               =gdkd$Gene %in% as.character(cnv[,1]) #Database genes mutated in patient
  druggable           =subset(gdkd,index)
  rownames(druggable) = NULL
  
  # Retain "copy number gain/loss" variants in the database
  index               =grep("copy number|promoter methylation or deletion",druggable[,"Description"])
  druggable           =druggable[sort(unique(index)),]   
  rownames(druggable) = NULL
  
  # Matching patient variants
  index=character()
  if (nrow(druggable) != 0){
    for (i in 1:nrow(druggable)){
      gene = as.character(druggable[i,"Gene"])
      pvar = as.character(cnv[which(cnv[,1]==gene),2])[1]
      
      # Get "mutation or copy number loss"
      if (grepl("copy number loss|promoter methylation or deletion",druggable[i,"Description"]) && pvar == "deletion"){
       druggable[i,"Patient_variant"] = pvar
      }
      if (grepl("copy number gain",druggable[i,"Description"]) && pvar == "amplification"){
        druggable[i,"Patient_variant"] = pvar
      }   
    }
  }
  index               = which(druggable$Patient_variant != "")
  druggable           = druggable[index,]          
  rownames(druggable) = NULL
  return(druggable)
}


match_TX_GDKD =function(tx,db = read.delim("data/GDKD.csv",sep="\t")){
  ## Given a list of gene rearrangements, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
  ## input: tx is a vector with gene rearrangements. Each gene rearrangement is defined by two dash-separated genes "gene1-gene2".
  ## output: returns those Gene Drug Knowledge Database rows matching to the input gene rearrangements.
  
  gdkd = db
  gdkd$Patient_variant = ""
  druggable                 = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
  colnames(druggable)       = colnames(gdkd)
  
  # Retain "fusion gene" or "rearrangement" variants
  index               = grep("fusion gene|rearrangement",gdkd[,"Description"])   
  druggable           = gdkd[index,]   
  rownames(druggable) = NULL

  # grep patient single gene to gdkd$Gene
  for (fusion in tx){
    genes  = strsplit(fusion,",|-")[[1]]
    for (gene in genes){
      i  = grep(gene,druggable$Gene,ignore.case=TRUE)
      druggable[i,"Patient_variant"] = fusion
    }
  }
  
  ## Retain rows with "Patient_variant"
  index=character()
  index= which(druggable$Patient_variant != "")
  druggable=druggable[index,]
  rownames(druggable) = NULL
  
  return(druggable)
 
}

match_WT_GDKD = function(snv, cnv, cancer_gdkd,db = read.delim("data/GDKD.csv",sep="\t")){
  ## Searches for wild type variants in Gene Drug Knowledge Database not matching any gene in the inputs SNVs or CNVs.
  ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## input: cancer_gdkd is a character vector with a cancer type used by GDKD.
  ## output: returns those GDKD wild type rows (annotated for disease=cancer_gdkd) not matching to the input SNVs and CNVs.

  gdkd = db
  gdkd$Patient_variant = ""
  druggable = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
  colnames(druggable)       = colnames(gdkd)
   
  wt  = grep("wild type",gdkd[,"Variant"],ignore.case=TRUE)
  wt  = wt[!(gdkd[wt,"Gene"] %in% c(as.character(cnv[,1]),as.character(snv[,1])))] # keep those that are not in cnv nor snv


  #check cancer type=cancer_gdkd
  same_cancer=c()
  for (ck in cancer_gdkd){
    same_cancer = c(which(gdkd[wt,"Disease"]==ck),same_cancer)
  }
  druggable  = gdkd[wt[unique(same_cancer)],]
  
  return(druggable)
  
}
################
## 2) CIVIC   ##
################

match_SNV_CIVIC = function(snv,db = read.delim("data/CIViC.csv",sep="\t")){
  ## Given a list of SNVs searches for matching variants in CIViC database.
  ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## output: returns those CIViC rows matching to the input SNVs.

  civic = db
  civic$Patient_variant = ""
  druggable             = data.frame(matrix(ncol=ncol(civic),nrow=0))
  colnames(druggable)   = colnames(civic)
  
  # Gene match
  index               = as.character(civic$gene) %in% as.character(snv[,1])  #Database genes mutated in patient
  druggable           = subset(civic,index)                                  #Subset from complete DB
  rownames(druggable) = NULL

 # Keep matching variants
  if (nrow(druggable) != 0){
    for (i in 1:nrow(druggable)){
      gene       = druggable[i,"gene"]
      p_var_type = as.character(snv[which(as.character(snv[,1])==gene),2])
      p_variants = unique(as.character(snv[which(as.character(snv[,1])==gene),3]))
      for (p_variant in p_variants){
        ### Any mutation
        if (druggable[i,"variant"]=="MUTATION" || druggable[i,"variant"]=="LOSS OF FUNCTION"){
          druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
        }
        ### frameshift,ins,del
        if (grepl("frameshift|frame shift",druggable[i,"variant"], ignore.case = TRUE) && grepl("frame_shift",p_var_type,ignore.case = TRUE) ){
          druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
        }
        ### stopgain ~ deletion
        if (grepl("deletion|loss",druggable[i,"variant"], ignore.case = TRUE) && grepl("nonsense",p_var_type,ignore.case = TRUE) ){
          druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=" "))
        }
        ### Missense & Specific variant
        if (grepl("missense",p_var_type[match(p_variant,p_variants)],ignore.case = TRUE)){
          if (grepl(p_variant,druggable[i,"variant"])){ 
            druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
          }
          ## Repurposing - New variants         
          else { 
            if (grepl(gene,"BRAF,KIT,FGFR3,KRAS,PDGFRA,EGFR,AKT1,MTOR,ALK",ignore.case = TRUE)){ 
              if (gene=="EGFR" &&  druggable[i,"variant"] != "T790M" && druggable[i,"variant"] != "C797S"){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
              if (gene=="MTOR" && druggable[i,"variant"] != "F2108L"){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
              if (gene=="AKT1" && druggable[i,"variant"] != "Q79K"){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
              if (gene=="BRAF" && grepl("V600",druggable[i,"variant"])==F){
                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
            }
            else {
              druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ") 
            }
          }
        }


      }
    }
  }

  ## Keep all those entries with "Patient_variant" and WT
  index      = which(druggable$Patient_variant != "")
  druggable  = druggable[index,]   #Subset from first selection
 # druggable  = rbind(druggable,civic[wt,])
                                
  return(druggable)
  
}


match_CNV_CIVIC =function(cnv,db = read.delim("data/CIViC.csv",sep="\t")){
  ## Given a list of CNVs searches for matching variants in CIViC database.
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## output: returns those CIViC rows matching to the input CNVs.

  civic = db
  civic$Patient_variant = ""
  druggable             = data.frame(matrix(ncol=ncol(civic),nrow=0))
  colnames(druggable)   = colnames(civic)

  # Gene match
  index     = as.character(civic$gene) %in% as.character(cnv[,1])  #Database genes mutated in patient
  druggable = subset(civic,index)                                  #Subset from complete DB
  rownames(druggable) = NULL

 # WT
 # wt  = grep("wild type",civic[,"variant"],ignore.case=TRUE)
 # wt  = wt[!(civic[wt,"gene"] %in% as.character(cnv[,1]))] # keep those that are not in cnv

  # Keep "amplification" or "deletion" variants
  index     = grep("amplification|^deletion$|loss",druggable[,"variant"],ignore.case=TRUE)
  druggable = druggable[sort(unique(index)),]     
  rownames(druggable) = NULL

  # Keep matching patient variants
  index=character()
  if (nrow(druggable) != 0){
    for (i in 1:nrow(druggable)){
      gene = as.character(druggable[i,"gene"])
      pvar = as.character(cnv[which(cnv[,1]==gene),2])[1]
      # Get "mutation or copy number loss"
      if (grepl("loss|deletion",druggable[i,"variant"],ignore.case=TRUE) && pvar == "deletion"){
        druggable[i,"Patient_variant"] = pvar
      }
      if (grepl("amplification",druggable[i,"variant"],ignore.case=TRUE) && pvar == "amplification"){
        druggable[i,"Patient_variant"] = pvar
      }
    }
  }
  index               = which(druggable$Patient_variant != "")
  druggable           = druggable[index,]         
 # druggable           = rbind(druggable,civic[wt,])   
  rownames(druggable) = NULL

  return(druggable)
}


match_WT_CIVIC  = function(snv,cnv,cancer_civic,db = read.delim("data/CIViC.csv",sep="\t")){
  ## Searches for wild type variants in CIViC not matching any gene in the inputs SNVs or CNVs.
  ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## input: cancer_civic is a character vector with a cancer type used by CIViC.
  ## output: returns those CIViC wild type rows (annotated for disease=cancer_civic) not matching to the input SNVs and CNVs.

  civic = db
  wt  = grep("wild type",civic[,"variant"],ignore.case=TRUE)
  wt  = wt[!(civic[wt,"gene"] %in% c(as.character(cnv[,1]),as.character(snv[,1])))] # keep those that are not in cnv nor snv

  #check cancer type=cancer_civic
  same_cancer=c()
  for (ck in cancer_civic){
    same_cancer = c(which(civic[wt,"disease"]==ck),same_cancer)
  }

  druggableCIVIC_wt  = civic[same_cancer,]
  
  return(druggableCIVIC_wt)
}



############################
## 3) TARGET              ##
## 4) MERIC-BERNSTAM DB   ##
############################

match_TARGET_MERIC = function(snv=c(),cnv=c(),tx=c(),db=read.delim("data/TARGET_MERIC.csv",sep="\t")){
  ## Given a list of SNVs, CNVs and gene rearrangements, searches for matching variants in TARGET database and gene list from Meric-Bernstam et al 2015.
   ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## input: tx is a vector with gene rearrangements. Each gene rearrangement is defined by two dash-separated genes "gene1-gene2".
  ## output: returns those TARGET + MERIC-BERNSTAM DB rows matching to the input SNVs, CNVs and gene rearrangements.

  target_meric = db
  druggable  = data.frame(matrix(ncol=5,nrow=0))
  druggable2 = data.frame(matrix(ncol=5,nrow=0))
  druggable3 = data.frame(matrix(ncol=5,nrow=0))
  druggable4 = data.frame(matrix(ncol=5,nrow=0))

  colnames(druggable)  = c(colnames(target_meric)[1:4],"Patient_variant")
  colnames(druggable2) = c(colnames(target_meric)[1:4],"Patient_variant")
  colnames(druggable3) = c(colnames(target_meric)[1:4],"Patient_variant")
  colnames(druggable4) = c(colnames(target_meric)[1:4],"Patient_variant")

  # Match SNVs
  if(is.null(snv) || nrow(snv)==0)warning('No SNVs provided!')
  else{
    index     = target_meric$Gene %in% as.character(snv[grep("missense",snv[,2],ignore.case=TRUE),1]) 
    druggable = subset(target_meric,index,1:4)
    druggable$Patient_variant = snv[na.omit(match(target_meric$Gene,as.character(snv[grep("missense",snv[,2],ignore.case=TRUE),1]))),3]
    index     = grepl("Mutation",druggable[,"Types_of_recurrent_alterations"],ignore.case = TRUE)
    druggable = subset(druggable,index)
  }

  # Match CNVs
  if(is.null(cnv) || nrow(cnv)==0)warning('No CNVs provided!')
  else{
    index      = target_meric$Gene %in% as.character(cnv[,1])
    druggable2 = subset(target_meric,index,1:4)
    druggable2$Patient_variant = cnv[na.omit(match(target_meric$Gene,as.character(cnv[,1]))),2]
    index      = grepl("amplification|deletion|biallelic",druggable2[,"Types_of_recurrent_alterations"],ignore.case = TRUE)
    druggable2 = subset(druggable2,index)
    # Match specific structural alteration
    if (nrow(druggable2)!=0){
      druggable2[,"Types_of_recurrent_alterations"] = gsub("[Ii]nactivation","deletion",druggable2[,"Types_of_recurrent_alterations"])
      index      = paste(lapply(1:nrow(druggable2), function(x) grep(druggable2[x,5],druggable2[x,3], ignore.case = TRUE)))
      index      = index==1
      druggable2 = druggable2[index,]
    }
  }
  # Match TX
  if(is.null(tx) || length(tx)==0)warning('No gene rearrangements provided!')
  else{    
    druggable3  = subset(target_meric,grepl("Rearrangement",target_meric[,"Types_of_recurrent_alterations"]),1:4)
    for (fusion in tx){
      genes  = strsplit(fusion,",|-")[[1]]
      for (gene in genes){
        i  = grep(gene,druggable3$Gene,ignore.case=TRUE)
        druggable3[i,"Patient_variant"] = paste(fusion,"fusion",sep=" ")
      }
    }
    index                = character()
    index                = which(druggable3$Patient_variant != "")
    druggable3           = druggable3[index,]  #Subset those with Patient variant
    rownames(druggable3) = NULL
  }

  # Match Nonsense mut - repurposing rule (nonsense ~ deletion, biallelic inactivation)
  if(!is.null(snv) && nrow(snv)!=0){
    index      = target_meric$Gene %in% as.character(snv[grep("nonsense",snv[,2],ignore.case=TRUE),1]) 
    druggable4 = subset(target_meric,index,1:4)
    druggable4$Patient_variant = snv[na.omit(match(target_meric$Gene,as.character(snv[grep("nonsense",snv[,2],ignore.case=TRUE),1]))),3]
    index      = grepl("deletion|inactivation",druggable4[,"Types_of_recurrent_alterations"],ignore.case = TRUE)
    druggable4 = subset(druggable4,index)
  }

  # rbind snv + cnv + tx druggable
  if(nrow(druggable)!=0 || nrow(druggable2)!=0 || nrow(druggable3)!=0 || nrow(druggable4)!=0){
    druggable = rbind(druggable,druggable2, druggable3, druggable4)
  }
  return(druggable[,c(1,5,2,3,4)]) 
}



