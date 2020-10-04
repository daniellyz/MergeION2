#' Substructure recommendation for unknown spectrum
#'
#' The function recommends substructures for an unknown spectrum 
#'
#' @param query_spectrum  Two-column data matrix. Two columns represent m/z and intensity of query tandem spectrum. At least 3 valid peaks should be provided
#' @param params.search General parameters for matching fragments to 
#' \itemize{
#'  \item{mz_search:}{ Numeric. Absolute mass tolerance in Da for fragment match.}
#'  \item{tops:}{ Integer. Top n substructures recommended.}
#' }
#' @param params.query.sp Parameters for matching an unknown spectra
#' \itemize{
#'   \item{prec_mz:}{ Numeric. Precursor mass of query spectrum. please put the correct value of precursor mass.}
#'   \item{use_prec:}{ Boolean. If set to TRUE, precursor mass is used to calculate the neutral loss.}
#'   \item{compound.type:}{ Character. Either "Metabolite" or "Drug". Different models are applied to Metabolite or Drug.}
#' }
#' 
#' @return Recommended substructures (ordered) as a two-column matrix: SMILES code and score
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#'
#' @examples
#' 
#' url_test = "https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/spectra/caffeine.txt"
#' 
#' test = read.table(url_test, header = F, sep=" ")
#' 
#' results = library_messar(query_spectrum = test, params.query.sp = list(prec_mz = 195.088, use_prec = T, compound.type = "Drug"))
#'
#' @export

library_messar<-function(query_spectrum =NULL, params.search = list(mz_search = 0.01, tops = 5),
                         params.query.sp = list(prec_mz = 100, use_prec = T, compound.type = "Metabolite")){

  options(stringsAsFactors = FALSE)
  options(warn=-1)
  output = NULL
  
  ####################
  ### Process inputs##
  ####################
  
  mz_search = params.search$mz_search
  ppm_search = params.search$ppm_search
  tops = params.search$tops
  
  prec_mz = as.numeric(params.query.sp$prec_mz)
  use_prec = params.query.sp$use_prec
  cpd_type = params.query.sp$compound.type
  
  if (!use_prec){prec_mz = 10000} # Not consider the precursor mass
  
  if (cpd_type=="Drug"){
    data(DRUG_RULE_DB)
    rules0 = drug_rules0
    rules_id_list = drug_rules_id_list
    rule_fragments = drug_rule_fragments
    rule_nloss = drug_rule_nloss
   # rules_training = drug_rules_training
  } else {data(METABOLITE_RULE_DB)}
  
  
  ###############################
  ### Preprocess query spectrum:#
  ###############################
  
  dat = query_spectrum[,1:2, drop = FALSE]
  dat = denoise_query_spectrum(dat, prec_mz, 200, 0.001)
  NP = nrow(dat)
  
  if (NP<3){
    message("The query spectrum must contain at least 3 valid peaks!")
    return(NULL)
  }
  
  frags = dat[,1]
  nloss = prec_mz - frags
  nloss = nloss[nloss>0.5] # Must be higher than 0.5
  
  ################
  ##### MESSAR####
  ################
  
  id_matched = c() # ID matched to the aligned features
  
  for (frag in frags){
    error = abs(rule_fragments$Mass-frag)
    valid = which.min(error)
    if (error[valid]<mz_search){
      id_matched = c(id_matched, rule_fragments$ID[valid])
    }
  }
  
  for (nl in nloss){
    error = abs(rule_nloss$Mass - nl)
    valid = which.min(error)
    if (error[valid]<mz_search){
      id_matched = c(id_matched, rule_nloss$ID[valid])
    }
  }

  # Match to rules and scoring:
    
  matched_rule_score = lapply(rules_id_list, function(x) intersect(as.character(id_matched), x))
  matched_rule_score = sapply(matched_rule_score,length)/log2(rules0$SIZE)

  # Scoring substructures:
    
  substructure_estimated = cbind.data.frame(SUBSTRUCTURE = rules0$SUBSTRUCTURE, SCORE = matched_rule_score)
  substructure_estimated = substructure_estimated[substructure_estimated[,2]>0,,drop=FALSE]
  substructure_estimated = substructure_estimated[order(substructure_estimated[,2], decreasing=T),]

  if (nrow(substructure_estimated)>0){
    tops1 = min(nrow(substructure_estimated), tops)
    output = substructure_estimated[1:tops1,]
  }
  
  return(output)
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
