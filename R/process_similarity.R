#' Searching a query spectrum in a spectral library
#'
#' The function calculates spectral similarity of a query spectrum to an existing spectral library
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#' 
#' @export

process_similarity<- function(query_spectrum, polarity = "Positive", prec_mz = 100, use.prec = FALSE, input_library = NULL,  
                              method = c("Precision", "Recall", "F1", "Cosine", "Spearman", "MassBank", "NIST"), 
                              prec_ppm_search = 10, frag_mz_search = 0.005, min_frag_match = 6){
  
  options(stringsAsFactors = FALSE)
  options(warn=-1)
  gc()

  if (!is.null(query_spectrum)){if (ncol(query_spectrum)<2){stop("Spectrum must have 2 columns m/z and intensity!")}}
  if (!(polarity %in% c("Positive", "Negative"))){stop("Polarity of query spectrum must be positive or negative")}
  if (min_frag_match<5){stop("min_frag_match should not be smaller than 5!")}

  ################################
  ### Preprocess query spectrum:##
  ################################
  
  prec_mz = as.numeric(prec_mz)
  
  dat = denoise_query_spectrum(query_spectrum, prec_mz, 500, 0.01)
  NP = nrow(dat)
  
  if (NP<3){
    message("The query spectrum must contain at least 3 valid peaks!")
    return(NULL)
  }
  
  ###########################################
  ### Read, filter and check input library###
  ###########################################
  
  input_library = library_reader(input_library)
  
  if (is.null(input_library$consensus)){stop("To allow spectral library search, consensus library must be created!")}
  
  db_profile = input_library$network$db_profile
  db_feature = input_library$network$db_feature
  consensus_library  = input_library$consensus
  
  if (use.prec){
    consensus_library = process_query(consensus_library, query = paste0("PEPMASS=", prec_mz), ppm_search = prec_ppm_search)
    consensus_library = consensus_library$SELECTED
  
    ID_SELECTED = consensus_library$metadata$ID
    valid1 = match(as.character(ID_SELECTED), colnames(db_profile))
    db_profile = db_profile[,valid1,drop=FALSE]

    valid2 = which(apply(db_profile, 1, function(x) sum(x>0))>0) # Remove all zero features
    db_profile = db_profile[valid2,,drop=FALSE]
    db_feature = db_feature[valid2,,drop=FALSE]
  }
  
  if (nrow(db_profile)==0){return(NULL)}
    
  ##############################################
  ### Pre-search common fragment/neutral loss###
  ##############################################
  
  db_profile1 = c()
  db_feature1 = c()
  dat1 = c()
  
  for (i in 1:NP){
    
    frags = dat[i,1]
    nls = prec_mz - dat[i,1]
    
    mze = abs(frags - db_feature[,2])
    nle = abs(nls - db_feature[,2])
  
    valid1 = which(mze <= frag_mz_search & db_feature$Type == "Frag")[1]
    valid2 = which(nle <= frag_mz_search & db_feature$Type == "Nloss")[1]
    
    valid = c(valid1, valid2)

    if (length(valid)>0){
      valid = valid[1]
      db_profile1 = rbind(db_profile1, db_profile[valid,,drop=FALSE])
      db_feature1 = rbind(db_feature1, db_feature[valid,,drop=FALSE])
      dat1 = rbind(dat1, dat[i,,drop=FALSE])
    }
  }
  
  if (is.null(db_profile1)){return(NULL)}
  
  if (!is.null(db_profile1)){
    if (nrow(db_profile1)==0){
    return(NULL)
  }}
  
  # Filter out db samples with fewer than minimum fragment matches:
  
  peak_matches = apply(db_profile1, 2, function(x) sum(x>0, na.rm = T))
  valid = which(peak_matches>=min_frag_match)  
  if (length(valid)==0){return(NULL)}
  db_profile1 = db_profile1[,valid,drop = FALSE]
  consensus_library1  = list(metadata = consensus_library$metadata[valid,,drop=FALSE], sp = consensus_library$sp[valid])
  
  # Filter out empty db features:
  
  feature_matches = apply(db_profile1, 1, function(x) sum(x>0, na.rm = T))
  valid = which(feature_matches>0)
  if (length(valid)==0){return(NULL)}
  
  db_profile = db_profile1[valid,,drop = FALSE]
  db_feature = db_feature1[valid,,drop = FALSE]
  dat = dat1[valid,,drop=FALSE]
  
  NDB = ncol(db_profile)

  ###########################
  ### Calculate Similarity###
  ###########################

  # Normalize first the spectra:
  
  dat[,2] = dat[,2]/max(dat[,2])*100 # Normalize
  db_profile <- apply(db_profile, 2, function(x) x/max(x)*100)

  # Calculate useful info:
  NP_query = nrow(dat) # Nb of peaks in query
  NP_reference = sapply(consensus_library1$sp, nrow)
  nb_matches = apply(db_profile, 2, function(x) sum(x>0))
  
  if (method == "Precision"){
    sim = nb_matches/NP_query
  }
    
  if (method == "Recall"){
    sim = nb_matches/NP_reference
  }
  
  if (method == "F1"){
    sim = 2*(nb_matches/NP_reference*nb_matches/NP_query)/(nb_matches/NP_reference+nb_matches/NP_query)
  }

  if (method == "Cosine"){
    sim = cor(dat[,2], db_profile, method = "pearson")/2 + 0.5
  }
  
  if (method == "Spearman"){
    sim = cor(dat[,2], db_profile, method = "spearman")/2 + 0.5
  }
  
  if (method == "MassBank"){
    dat_weighted = (db_feature$Mass^2)*(dat[,2]^0.5)
    db_profile_weighted = (db_feature$Mass^2)*(db_profile^0.5)
    sim = cor(dat_weighted, db_profile_weighted, method = "pearson")/2 + 0.5
  }
  
  if (method == "NIST"){
    dat_weighted = db_feature$Mass*dat[,2]
    db_profile_weighted = db_feature$Mass*db_profile
    sim = cor(dat_weighted, db_profile_weighted, method = "pearson")/2 + 0.5
  }
  
  sim = round(as.numeric(sim),2)
  sim.scores = cbind.data.frame(ID = colnames(db_profile), PEAK.MATCHES = nb_matches, SCORES = sim)

  ########################################
  ### Filter, add formula and output #####
  ########################################
  
  sim.scores = sim.scores[sim.scores$SCORES>0,,drop=FALSE]
  if (nrow(sim.scores)==0){return(NULL)}
  
  if ("FORMULA" %in% colnames(consensus_library$metadata)){
    all_formulas = consensus_library$metadata$FORMULA
    valid = match(sim.scores[,1], consensus_library$metadata$ID)
    sim.scores$FORMULA = all_formulas[valid]
  } else {sim.scores$FORMULA = "N/A"}
  
  sim.scores$ID[sim.scores$ID==""] = "N/A"
  sim.scores$FORMULA[sim.scores$FORMULA==""] = "N/A"
  sim.scores = sim.scores[order(sim.scores$SCORES, decreasing = T),]
  
  return(sim.scores)
}

##########################
### Internal functions:###
##########################

denoise_query_spectrum<-function(sp, mz0, max_peak, min_relative){
  
  denoised_spectrum = matrix(c(0,0),1,2)
  
  if (nrow(sp)>0){

    # Filter top peaks:
    
    sp = sp[order(sp[,2], decreasing = T),,drop=FALSE]
    tops = min(max_peak, nrow(sp))  
    sp = sp[1:tops,,drop=FALSE]
    
    # Normalize to 100:
    
    sp1 = sp
    sp1[,2] = sp1[,2]/max(sp1[,2])*100
    
    # Relative Intensity filter:
    
    filter = which(sp1[,2]>=min_relative & sp1[,1]<mz0-1)
    sp = sp1
    sp = sp[filter,,drop=FALSE]
    
    # Check validity:
    
    if (nrow(sp)>0){
      sp = sp[order(sp[,1]),,drop=FALSE]
      denoised_spectrum = sp
    }
  }
  return(denoised_spectrum)
}