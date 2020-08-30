#' Create pseudo-library file from a LC-MS/MS data file through XC-MS centWave non-targeted feature screening
#'
#' Function used by library_generator to screen LC-MS features from unknown mzML/mzXML file
#'
#' @param raw_data_file Character. The LC-MS/MS file name from which fragmented LC-MS features are extracted. All files must have be in centroid-mode with mzML or mzMXL extension!
#' @param polarity character. Either "Positive" or "Negative". Ion mode of the LC-MS/MS file. 
#' @param ppm_search Numeric. Absolute mass tolerance in ppm for feature screening and extraction
#' @param rt_search Numeric. Absolute retention tolerance in seconds for feature screening and extraction
#' @param baseline Numeric. Absolute intensity threshold above which is considered as a mass peak, otherwise noise.
#' 
#' @importFrom xcms xcmsSet
#' @importFrom CAMERA xsAnnotate groupFWHM findIsotopesWithValidation getPeaklist
#' @importFrom tools file_path_sans_ext
#' @importFrom MSnbase readMSData rtime precursorMz
#' @importFrom proxy dist
#' @importFrom stats hclust cutree
#' 
#' @export
#'
#'
process_dda<-function(raw_data_file, polarity = c("Positive", "Negative"),
                      ppm_search = 10, rt_search =12, baseline = 1000){
  
  options(stringsAsFactors = F)
  options(warn=-1)
  
  #################
  ### Initialize ##
  #################
  
  if (!is.character(raw_data_file)){
    stop("Please provide a chromatogram file!")}
  
  if (!(file_ext(raw_data_file) %in% c("mzML", "mzXML"))){
    stop("The chromatogram file must be in mzML or mzXML format!")}
  
  if (!(polarity %in% c("Positive", "Negative"))){
    stop("Polarity must be either Positive or Negative")
  }
  
  feature_mz = c()
  feature_rt = c()
  feature_type = c()
  
  ###################
  ### MS2 Screener ##
  ###################
  
  MS2_Janssen <- try(readMSData(raw_data_file, msLevel = 2,  verbose = FALSE),silent=T)
  
  ### Extract feature m/z
  
  feature_mz = c(feature_mz, precursorMz(MS2_Janssen)) # Label precursor mass
  feature_rt = c(feature_rt, rtime(MS2_Janssen)/60) # RT in min
  feature_type = c(feature_type, rep("DDA", length(precursorMz(MS2_Janssen))))
  
  ##############################################
  ### Extract features from MS1 scans with XCMS#
  ##############################################
  
  xs <- xcmsSet(raw_data_file, method = "centWave", ppm = ppm_search*2, noise = baseline,
                snthresh = 30, prefilter = c(6, baseline))
  
  ### Remove isotopes
  
  an <- xsAnnotate(xs)
  an <- groupFWHM(an, perfwhm=2)
  an <- findIsotopesWithValidation(an, ppm = ppm_search,  mzabs = 0.01, maxcharge = 2, database = "pubchem")  # optional but recommended.
  
  peaks = getPeaklist(an)
  peaks = data.frame(peaks)
  
  iso1 = grep("[M+1]+", peaks$isotopes, fixed = T)
  iso2 = grep("[M+2]+", peaks$isotopes, fixed = T)
  peaks = peaks[-c(iso1, iso2),]
  
  ### Extract feature m/z
  
  feature_mz = c(feature_mz, peaks$mz)
  feature_rt = c(feature_rt, peaks$rt/60) # RT in min
  feature_type = c(feature_type, rep("XCMS", length(peaks$mz)))
  
  #####################################
  ### Aggregate extracted m/z and RT ##
  #####################################
  
  prec_dist = proxy::dist(feature_mz, ppm_distance1)
  prec_clust = hclust(prec_dist, method = "average")
  prec_tree = cutree(prec_clust, h = ppm_search*2)
  
  rt_dist = proxy::dist(feature_rt, rt_distance1)
  rt_clust = hclust(prec_dist, method = "average")
  rt_tree = cutree(prec_clust, h = rt_search*2/60)
  
  prec_rt_label = paste(prec_tree, rt_tree, sep = "-")
  prec_rt_class = unique(prec_rt_label)
  
  MS2_PEPMASS = c()
  MS2_RT = c()
  
  for (pt in prec_rt_class){
    
    temp_ind = which(prec_rt_label == pt)
    temp_mz = median(feature_mz[temp_ind])
    temp_rt = median(feature_rt[temp_ind])
    temp_type = feature_type[temp_ind]
    
    if (sum(temp_type=="XCMS")>0 & sum(temp_type=="DDA")>0){
      MS2_PEPMASS = c(MS2_PEPMASS, temp_mz)
      MS2_RT = c(MS2_RT, temp_rt)
    }
  }
  
  #######################
  ### Create metadata ###
  #######################
  
  MS2_metadata = matrix("N/A", length(MS2_PEPMASS), 5)
  colnames(MS2_metadata) = c("ID", "PEPMASS","RT", "ADDUCT", "FILENAME")
  MS2_metadata = data.frame(MS2_metadata)
  
  MS2_metadata$PEPMASS = MS2_PEPMASS
  MS2_metadata$RT = MS2_RT
  MS2_metadata$ID = paste0(basename(file_path_sans_ext(raw_data_file)),"_",1:nrow(MS2_metadata))
  MS2_metadata$FILENAME = basename(raw_data_file)
  
  if (polarity=="Positive"){MS2_metadata$ADDUCT = "M+H"}
  if (polarity=="Negative"){MS2_metadata$ADDUCT = "M-H"}
  
  return(MS2_metadata)
 
}

####################################
###########Internal function########
####################################

ppm_distance1<-function(x,y){
  x = as.numeric(x)
  y = as.numeric(y)
  ppm = abs((x-y)/y*1000000)
  if (x<120 & abs(x-y)<0.01){
    ppm = 1
  }
  return(ppm)
}

rt_distance1<-function(x,y){
  return(abs(x-y))
}