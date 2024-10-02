#' Create additional metadata for library generation
#'
#' Function used by library_generator() to improve and adapt input metadata
#' 
#' @importFrom rcdk parse.smiles get.mol2formula
#' @importFrom OrgMassSpecR MonoisotopicMass ListFormula
#' @importFrom webchem is.smiles
#' @importFrom tools file_ext
#' @importFrom stringr str_detect
#' @inheritParams library_generator
#' 
#' @export
#
process_metadata<-function(ref, processing.algorithm = c("Default", "compMS2Miner", "RMassBank"),
                          polarity = c("Positive", "Negative"), add.adduct = T, adductType = NULL){
  
  ref0 = ref # Backup

  new_metadata = c()
  ref_adducts = NULL
    
	 if(is.null(adductType)){
      if (polarity == "Positive"){ref_adducts = c("M+H","M+2H","M+Na","M+K","M+NH4", "M+", "M+H-H20")}
      if (polarity == "Negative"){ref_adducts = c("M-H","M+Cl", "M+HCOO-","M+CH3COO-")}
  } else {
	    ref_adducts <- adductType
  }

  ##########################
  ### Complete metadata ####
  ##########################
  
  ref$ID=as.character(ref$ID) # Make sure IDs are characters
  
  # Replace invalid pepmass and retention by N/A
  
  if ("PEPMASS" %in% colnames(ref)){
    NR = 1:nrow(ref)
    ref_mz = as.numeric(ref$PEPMASS)
    nonvalid = setdiff(NR, which(ref_mz>50 & ref_mz<3000))
    ref$PEPMASS[nonvalid] = "N/A"
  } else {ref$PEPMASS = "N/A"}
  
  if ("RT" %in% colnames(ref)){
    NR = 1:nrow(ref)
    ref_rt = as.numeric(ref$RT)
    nonvalid = setdiff(NR, which(ref_rt>0 & ref_rt<60))
    ref$RT[nonvalid] = "N/A"
  } else {ref$RT = "N/A"}
  
  # Remove opposite polarity and add default polarity to metadata:
  
  if ("IONMODE" %in% colnames(ref)){
    if (polarity=="Positive"){valid = which(ref$IONMODE=="Positive")}
    if (polarity=="Negative"){valid = which(ref$IONMODE=="Negative")}
    ref = ref[valid,,drop=FALSE]
  }
  
  ref$IONMODE = polarity
  
  # Remove invalid adduct and add default M+H or M-H:
  
  if ("ADDUCT" %in% colnames(ref)){
    valid = which(ref$ADDUCT %in% ref_adducts)
    ref = ref[valid,,drop=FALSE]
  } else { # Otherwise all considered as default
    if (polarity=="Positive"){ref$ADDUCT = "M+H"}
    if (polarity=="Negative"){ref$ADDUCT = "M-H"}
  }

  # Set default charge states:
  
  ref$CHARGE = 1
  ref$CHARGE[ref$ADDUCT=="M+2H"] = 2
  
  # Check smiles code if provided:
  
  if ("SMILES" %in% colnames(ref)){
    nonvalid = which(!sapply(ref$SMILES,is.smiles))
    ref$SMILES[nonvalid] = "N/A"
  } else {ref$SMILES = "N/A"}
  
  if (!("FORMULA" %in% colnames(ref))){ref$FORMULA = "N/A"}
  
  # Add CAS and Names for RMassBank:
  
  if (processing.algorithm=="RMassBank"){
    ref$CAS = "N/A"
    ref$Name = "N/A"
    
    #for (i in 1:nrow(ref)){
    #kkk = cas_name_generator(ref$SMILES[i])
    #ref$CAS[i] = kkk$CAS
    #ref$Name[i] = kkk$Name
    #}
  }
  
  # Check filename:
  
  if ("FILENAME" %in% colnames(ref)){
    nonvalid = which(!(toupper(file_ext(ref$FILENAME)) %in% c("MZML", "MZXML", "CDF")))
    ref$FILENAME[nonvalid] = "N/A"
  } else {ref$FILENAME = "N/A"}
  
  # Calculate PEPMASS from smiles if PEPMASS not available:
  
  for (i in 1:nrow(ref)){
    k = NULL
    if (ref$SMILES[i]!="N/A" & ref$PEPMASS[i]=="N/A"){
      k = metadata_from_smiles(ref$SMILES[i], ref$ADDUCT[i])
    }
    if (!is.null(k)){
      ref$PEPMASS[i] = k$pepmass
      ref$FORMULA[i] = k$formula
      ref$THEOMASS[i] = k$theomass
    }
  }
  
  # Expand adduct if required by users:
  
  ref1 = ref[which(!ref$ADDUCT %in% c("M-H", "M+H")),,drop=FALSE]
  ref2 = ref[which(ref$ADDUCT %in% c("M-H", "M+H")),,drop=FALSE]
  
  if(nrow(ref2) != 0){
    ref2_new = expand_ref_adduct(ref2, polarity, add.adduct, adapt.smiles = processing.algorithm == "RMassBank")
  }else{
	  ref2_new = NULL
  }
  
  ref = rbind.data.frame(ref1, ref2_new)
  
  # Change column order:
  
  ref1 = ref[, c("PEPMASS", "RT", "IONMODE", "ADDUCT", "CHARGE", "ID", "SMILES", "FORMULA", "FILENAME")]
  ref2 = ref[,which(!(colnames(ref) %in% c("PEPMASS", "RT", "IONMODE", "ADDUCT", "CHARGE", "ID", "SMILES", "FORMULA", "FILENAME"))), drop=FALSE]
  ref = cbind.data.frame(ref1, ref2)
  
  ###########################################
  ## Filter metadata according to algorithm##
  ###########################################
  
  if (processing.algorithm=="Default"){
    valid = which(ref$PEPMASS!="N/A")
    ref$PEPMASS = as.numeric(ref$PEPMASS)
  }
  
  if (processing.algorithm=="compMS2Miner"){
    valid = which(ref$PEPMASS!="N/A" & ref$RT!="N/A")
    ref$PEPMASS = as.numeric(ref$PEPMASS)
    ref$RT = as.numeric(ref$RT)
  }
  
  if (processing.algorithm=="RMassBank"){
    valid = which(!str_detect(ref$SMILES, "N/A") & ref$FILENAME!="N/A")
  }
  
  #######################
  # Combine final result#
  #######################
  
  if (length(valid)>0){new_metadata = ref[valid,,drop=FALSE]}
  
  return(new_metadata)
}

##########################
####1:Structure mining####
##########################

cas_name_generator<-function(smiles){
  
  options(stringsAsFactors = F)
  options(warn = -1)
  Name = "N/A"
  CAS = "N/A"
  
  base_url = "https://cactus.nci.nih.gov/chemical/structure/"
  url1 = paste0(base_url, smiles, "/names")
  url2 = paste0(base_url, smiles, "/cas")
  temp1 = try(as.character(readLines(url1)[1]), silent = T)
  temp2 = try(as.character(readLines(url2)[1]), silent = T)
  
  if (class(temp1)!="try-error"){Name = temp1}
  if (class(temp2)!="try-error"){CAS = temp2}
  return(list(Name = Name, CAS = CAS))
}

metadata_from_smiles<-function(smiles, adduct){
  
  output = NULL
  smiles1  = smiles
  
  mol1 = try(parse.smiles(smiles1)[[1]], silent = T)
  formula1 = try(get.mol2formula(mol1, charge = 0)@string, silent = T)
  list1 = try(ListFormula(formula1), silent = T)
  
  if (class(formula1) != "try-error"){
    
    NM = MonoisotopicMass(formula=list1) 
    
    if (adduct=="M+H"){pepmass = NM +1.007276}
    if (adduct=="M+2H"){pepmass = (NM + 1.007276*2)/2}
    if (adduct=="M+Na"){pepmass = NM +22.989221}
    if (adduct=="M+K"){pepmass = NM +38.963158}
    if (adduct=="M+NH4"){pepmass = NM +18.033826}
    if (adduct=="M+"){pepmass = NM +1.007276 - 1.007825}
    if (adduct=="M+H-H2O"){pepmass = NM -17.0033} 
    
    if (adduct=="M-H"){pepmass = NM -1.007276}
    if (adduct=="M-2H"){pepmass = (NM - 1.007276*2)/2}
    if (adduct=="M+Cl"){pepmass = NM +34.969401}
    if (adduct=="M+CH3COO-"){pepmass = NM + 59.013853}
    if (adduct=="M+HCOO-"){pepmass = NM + 44.998203}
    if (adduct=="M-"){pepmass = NM -1.007276+ 1.007825}
    
    output = list(smile = smiles1, formula = formula1, pepmass = pepmass, theomass = NM)
  }
  
  return(output)  
}

expand_ref_adduct <-function(ref, polarity, add.adduct = T, adapt.smiles = T){
  
  # The function expand the metadata by adding adducts
  
  if (polarity == "Positive"){
    
    NM = as.numeric(ref$PEPMASS) - 1.007276
    
    ref_plus_2H = ref
    ref_plus_2H$PEPMASS = (NM + 1.007276*2)/2  
    ref_plus_2H$ADDUCT = "M+2H"
    ref_plus_2H$CHARGE = 2
    if (adapt.smiles){ref_plus_2H$SMILES  = paste0(ref_plus_2H$SMILES,".[H].[H]")}
    
    ref_plus_na = ref
    ref_plus_na$PEPMASS = NM + 22.989221 
    ref_plus_na$ADDUCT = "M+Na"
    ref_plus_na$CHARGE = 1
    if (adapt.smiles){ref_plus_na$SMILES  = paste0(ref_plus_na$SMILES,".[Na]")}
    
    ref_plus_k = ref
    ref_plus_k$PEPMASS = NM + 38.963158 
    ref_plus_k$ADDUCT = "M+K"
    ref_plus_k$CHARGE = 1
    if (adapt.smiles){ref_plus_k$SMILES  = paste0(ref_plus_k$SMILES,".[K]")}
    
    ref_plus_nh4 = ref
    ref_plus_nh4$PEPMASS = NM + 18.033826 
    ref_plus_nh4$ADDUCT = "M+NH4"
    ref_plus_nh4$CHARGE = 1
    if (adapt.smiles){ref_plus_nh4$SMILES  = paste0(ref_plus_nh4$SMILES,".[NH4]")}
    
    ref_plus = ref
    ref_plus$PEPMASS = NM + 1.007276 - 1.007825 
    ref_plus$ADDUCT = "M+"
    ref_plus$CHARGE = 1
    
    if (adapt.smiles){ref$SMILES  = paste0(ref$SMILES,".[H]")}
    
    if (add.adduct){
      ref = rbind.data.frame(ref, ref_plus_2H, ref_plus_na, ref_plus_k, ref_plus_nh4, ref_plus)
    }}
  
  if (polarity == "Negative"){
    
    NM = as.numeric(ref$PEPMASS) + 1.007276
    
    ref_plus_cl = ref
    ref_plus_cl$PEPMASS = NM + 34.969401 
    ref_plus_cl$ADDUCT = "M+Cl"
    ref_plus_cl$CHARGE = 1
    if (adapt.smiles){ref_plus_cl$SMILES  = paste0(ref_plus_cl$SMILES,".[HCl]")}
    
    ref_plus_HCOO = ref
    ref_plus_HCOO$PEPMASS = NM + 44.998203 
    ref_plus_HCOO$ADDUCT = "M+HCOO-"
    ref_plus_HCOO$CHARGE = 1
    if (adapt.smiles){ref_plus_HCOO$SMILES  = paste0(ref_plus_HCOO$SMILES,".[H2CO2]")}
    
    ref_plus_CH3COO = ref
    ref_plus_CH3COO$PEPMASS = NM + 59.013853 
    ref_plus_CH3COO$ADDUCT = "M+CH3COO-"
    ref_plus_CH3COO$CHARGE = 1
    if (adapt.smiles){ref_plus_CH3COO$SMILES  = paste0(ref_plus_CH3COO$SMILES,".[C2H4O2]")}
    
    if (add.adduct){
      ref = rbind.data.frame(ref, ref_plus_cl, ref_plus_HCOO, ref_plus_CH3COO)
    }}
  
  return(ref)
}
