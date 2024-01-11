#' Searching a query spectrum in a spectral library
#' @param query_spectrum two columns spectrum list
#' @param use.loss boolean if natrual loss also be used, \code{prec_mz} can not be missing if set true.
#' @inheritParams library_query
#' The function calculates spectral similarity of a query spectrum to an existing spectral library
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#' 
#' @export

process_similarity<- function(query_spectrum, polarity = "Positive", prec_mz = 100, use.loss = TRUE, use.prec = FALSE, input_library = NULL,  
                              method = c("Precision", "Recall", "F1", "Cosine", "Spearman", "MassBank", "NIST", "HM", "Entropy"), 
                              prec_ppm_search = 10, frag_mz_search = 0.005, min_frag_match = 6){
  
  options(stringsAsFactors = FALSE)
  options(warn=-1)
  gc()

  if (!is.null(query_spectrum)){if (ncol(query_spectrum)<2){stop("Spectrum must have 2 columns m/z and intensity!")}}
  if (!(polarity %in% c("Positive", "Negative"))){stop("Polarity of query spectrum must be positive or negative")}
  if (min_frag_match<2){stop("min_frag_match should not be smaller than 2!")}

  ################################
  ### Preprocess query spectrum:##
  ################################
  
  prec_mz = as.numeric(prec_mz)
  if(is.na(prec_mz)) use.loss <- FALSE

  
  dat = denoise_query_spectrum(query_spectrum, prec_mz, 500, 0.01)
  NP = nrow(dat)
  dat0 = dat
  
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
  
  if (use.prec & !is.na(prec_mz)){
    
    consensus_library = process_query(consensus_library, query = paste0("PEPMASS=", prec_mz), ppm_search = prec_ppm_search)
    consensus_library = consensus_library$SELECTED

    ID_SELECTED = consensus_library$metadata$ID
    valid1 = match(as.character(ID_SELECTED), colnames(db_profile))
    db_profile = db_profile[,valid1,drop=FALSE]

    valid2 = which(apply(db_profile, 1, function(x) sum(x>0))>0) # Remove all zero features
    db_profile = db_profile[valid2,,drop=FALSE]
    db_feature = db_feature[valid2,,drop=FALSE]
  }
  
  if (nrow(db_profile)==0){
    return(NULL)
  }
  
  #else {db_profile <- apply(db_profile, 2, function(x) x/max(x)*100)}
    
  ##############################################
  ### Pre-search common fragment/neutral loss###
  ##############################################
  
  db_profile1 = c()
  db_feature1 = c()
  dat1 = c()
  
  for (i in 1:NP){

    frags = dat[i,1, drop = TRUE]
    mze = abs(frags - db_feature[,2])

    valid1 = intersect(which(mze <= frag_mz_search), which(db_feature$Type == "Frag"))
   
	if(use.loss){
		nls = prec_mz - dat[i,1, drop = TRUE]
		nle = abs(nls - db_feature[,2])
		
    valid2 = intersect(which(nle <= frag_mz_search), which(db_feature$Type == "Nloss"))  
    valid = c(valid1, valid2)
}else{
	valid = valid1
}
    if (length(valid)>0){
      tmp_profile =  apply(db_profile[valid,,drop=FALSE], 2, max)
      for (rk in 1:length(valid)){
        dat1 = rbind(dat1, dat[i,,drop=FALSE])
        db_feature1 = rbind(db_feature1, db_feature[valid[rk],,drop = FALSE])
        db_profile1 = rbind(db_profile1, tmp_profile)
      }
    }
  }
  
  rownames(db_profile1) = rownames(db_feature1)
  
  if (is.null(db_profile1)){return(NULL)}
  
  if (!is.null(db_profile1)){
    if (nrow(db_profile1)==0){
    return(NULL)
  }}
  
  # Filter out duplicated features
  
  filter = which(!duplicated(dat1[,1]))
  dat1 = dat1[filter,,drop=FALSE]
  db_profile1 = db_profile1[filter,,drop=FALSE]
  db_feature1 = db_feature1[filter,,drop=FALSE]
  
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

  # Calculate useful info:
  
  NP_query = nrow(dat) # Nb of peaks in query
  NP_reference = sapply(consensus_library1$sp, nrow)
  NP_reference = as.numeric(as.character(NP_reference)) # Convert NULL to NA

  nb_matches = apply(db_profile, 2, function(x) sum(x>0))
  nb_matches = as.numeric(as.character(nb_matches))

  if (method == "Precision"){
    sim = nb_matches/NP_query
  }
    
  if (method == "Recall"){
    sim = nb_matches/NP_reference
  }
  
  if (method == "F1"){
    sim = 2*(nb_matches/NP_reference*nb_matches/NP_query)/(nb_matches/NP_reference+nb_matches/NP_query)
  }

  if (method == "HM"){
    # Normalize first the spectra to sum!
    dat[,2] = dat[,2]/sum(dat[,2]) # Normalize
    db_profile <- apply(db_profile, 2, function(x) x/sum(x))
    sim = 2*colSums((dat[,2]*db_profile)/(dat[,2]+db_profile))
  }
  
  if (method == "Dot"){
    dat[,2] = dat[,2]/sum(dat[,2]) # Normalize
    db_profile <- apply(db_profile, 2, function(x) x/sum(x))
    sim = colSums((dat[,2]*db_profile)^2)/(sum((dat[,2]^2))*colSums(db_profile^2))
    sim = sim/2+0.5  
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
  
  if (method == "Entropy"){ 
    dat[,2] = dat[,2]/sum(dat[,2]) # Normalize
    
    global_min = min(apply(db_profile, 2, function(myvector) min(myvector[myvector > 0])))
    db_profile_imputed = db_profile
    db_profile_imputed[db_profile_imputed==0] = global_min
    db_profile_imputed <- apply(db_profile_imputed, 2, function(x) x/sum(x))
  
    mixed = (dat[,2]+db_profile_imputed)/2 # 1:1 Mixed
    mixed = apply(mixed, 2, function(x) x/sum(x)) # Normalize mixed again
    
    S_A = -sum(dat[,2]*log(dat[,2], base = exp(1)))
    #if (S_A<3){S_A = -sum(dat[,2]^(0.25+S_A*0.25)*log(dat[,2]^(0.25+S_A*0.25), base = exp(1)))}
    
    S_B = -colSums(db_profile_imputed*log(db_profile_imputed, base = exp(1)))
    #I_B_bis = db_profile_imputed[,S_B<3,drop=FALSE]^(0.25+S_B[S_B<3]*0.25)
    #S_B[S_B<3] =  -colSums(I_B_bis*log(I_B_bis, base = exp(1)))
    
    S_AB = -colSums(mixed*log(mixed, base = exp(1)))
    #I_AB_bis = mixed[,S_AB<3,drop=FALSE]^(0.25+S_AB[S_AB<3]*0.25)
    #S_AB[S_AB<3] = -colSums(I_AB_bis*log(I_AB_bis, base = exp(1)))
    
    sim  = 1 - (2*S_AB - S_A - S_B)/log(4, base = exp(1))
  }
  
  #if (is.na(sim)){sim = 0}
  #if (sim>1){sim = 1}
  
  sim[is.na(sim)] = 0
  sim[sim>1] = 1
  
  sim = round(as.numeric(sim),2)
  sim.scores = cbind.data.frame(ID = colnames(db_profile), PEAK.MATCHES = nb_matches, SCORES = sim)

  ########################################
  ### Filter, add formula and output #####
  ########################################
  
  sim.scores = sim.scores[sim.scores$SCORES>=0,,drop=FALSE]
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
    
		filter =if(is.na(mz0)) which(sp1[,2]>=min_relative) else which(sp1[,2]>=min_relative & sp1[,1]<mz0-1)
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
