#' Substructure recommendation for unknown spectrum
#'
#' The function recommends substructures for an unknown spectrum 
#'
#' @param query_spectrum  Two-column data matrix. Two columns represent m/z and intensity of query tandem spectrum. At least 3 valid peaks should be provided
#' @param prec_mz Numeric. Precursor mass of query spectrum (strongly recommended). Please set to 0 if unknown.
#' @param type Character. Two algorithms focused either on "drug" or "metabolite"
#' @param tops Integer. Top recommended substructures
#' @param max_peaks Integer. Top most intense peaks kept in the spectra
#' @param min_relative Numeric between 0.1 and 1. Minimum relative intensity of peaks to not be considered as noise.
#'
#' @return Recommended substructures (ordered) as a two-column matrix: SMILES code and score
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#'
#' @examples
#' 
#' url_test = "https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/caffeine.txt"
#' 
#' test = read.table(url_test, header = F, sep='\t')
#' 
#' results = library_messar_simple(query_spectrum = test, prec_mz = 195.088, type = "drug", tops = 3)
#'
#' @importFrom MSnbase fData readMgfData
#' @importFrom tools file_ext
#' @importFrom stringr str_replace_all fixed
#'
#' @export

library_messar_simple<-function(query_spectrum =NULL, prec_mz = 0, type = c("drug", "metabolite"), 
                                tops = 5, max_peaks = 200, min_relative = 0.1){
  
  options(stringsAsFactors = FALSE)
  options(warn=-1)
  output = NULL
  
  if (type=="drug"){
    data(DRUG_RULE_DB)
    rules0 = drug_rules0
    rules_id_list = drug_rules_id_list
    rule_fragments = drug_rule_fragments
    rule_nloss = drug_rule_nloss
    rules_training = drug_rules_training
  } else {data(METABOLITE_RULE_DB)}
  
  #################
  ### Check inputs:
  #################
  
  if (is.null(query_spectrum)){
    stop("Please provide a 2 column query spectrum!")
  } else {
    if (ncol(query_spectrum)<2){
      stop("Query spectrum must have at least 2 columns: m/z and intensity!")
    }
    if (nrow(query_spectrum)<3){
      stop("Query spectrum must have at least 3 peaks!")
    }
  }
  
  ###############################
  ### Preprocess query spectrum:
  ###############################
  
  dat = query_spectrum[,1:2, drop = FALSE]
  XP = nrow(dat)
  
  # Top peaks, normalize, cut only masses smaller than precursor and filter background noise:
  
  dat = dat[order(query_spectrum[,2], decreasing = T),] # Filter>0.1
  dat = dat[1:min(max_peaks, XP),]
  dat[,2] = dat[,2]/max(dat[,2])*100
  dat = dat[dat[,2]>min_relative,,drop = FALSE]
  dat = dat[order(dat[,1]),] # Filter>0.1
  
  if (prec_mz==0){prec_mz = 10000} # Not consider the precursor mass
  
  dat = dat[dat[,1]<prec_mz - 0.5,]
  
  if (nrow(dat)<3){
    stop("The query spectrum must contain at least 3 valid peaks!")
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
    if (error[valid]<0.003){
      id_matched = c(id_matched, rule_fragments$ID[valid])
    }
  }
  
  for (nl in nloss){
    error = abs(rule_nloss$Mass - nl)
    valid = which.min(error)
    if (error[valid]<0.003){
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
 
