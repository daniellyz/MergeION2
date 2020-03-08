#' Create pseudo-metadata for library generation
#'
#' Function used by library_generator to extract LC-MS features from MS1 chromatogram of DDA/targeted-mode
#'
#' @importFrom xcms xcmsSet
#' @importFrom CAMERA xsAnnotate groupFWHM findIsotopesWithValidation getPeaklist
#' @importFrom tools file_path_sans_ext
#' @importFrom MSnbase readMSData rtime tic fData readMgfData precursorMz polarity chromatogram
#' @importFrom proxy dist
#' @importFrom stats hclust cutree
#' @importFrom pracma trapz
#' @export
#'
metadata_MS1_screener<-function(raw_data_file, screener = c("XCMS", "DDA"), ref = NULL, max.charge = 1,
                                ppm_search = 20, rt_search =12, baseline = 1000, snthreshold = 30, min.points = 8){

  options(stringsAsFactors = F)
  options(warn=-1)

  ##############################################
  ### Extract features from MS1 scans with XCMS
  ##############################################

  if (screener == "XCMS"){

    xs <- xcmsSet(raw_data_file, method = "centWave", ppm = ppm_search*2,
                snthresh = snthreshold, prefilter = c(min.points, baseline))

    ### Remove isotopes

    an <- xsAnnotate(xs)
    an <- groupFWHM(an, perfwhm=2)
    an <- findIsotopesWithValidation(an, ppm = ppm_search,  mzabs = 0.01, maxcharge = max.charge, database = "pubchem")  # optional but recommended.

    peaks = getPeaklist(an)
    peaks = data.frame(peaks)
    
    iso1 = grep("[M+1]+", peaks$isotopes, fixed = T)
    iso2 = grep("[M+2]+", peaks$isotopes, fixed = T)
    iso3 = grep("[M+3]+", peaks$isotopes, fixed = T)
    iso4 = grep("[M+4]+", peaks$isotopes, fixed = T)
    peaks = peaks[-c(iso1, iso2, iso3, iso4),]

    if (xs@phenoData$class=="Positive"){
      polarity = "Positive"
      adduct = "M+H"
    } 
    
    if (xs@phenoData$class=="Negative"){
      polarity = "Negative"
      adduct = "M-H"
    }
    
    ID_list = paste0(basename(file_path_sans_ext(raw_data_file)),"_",1:nrow(peaks))
    MS1_metadata = cbind(peaks$mz, peaks$rt/60, polarity, adduct, 1, ID_list, peaks$intb)

    MS1_metadata = data.frame(MS1_metadata)
    colnames(MS1_metadata) = c("PEPMASS","RT","IONMODE","ADDUCT","CHARGE","ID","MS1_ABUNDANCE")
    tmp = c(1,2,5,7)
    MS1_metadata[,tmp] <- apply(MS1_metadata[,tmp],2,function(x) as.numeric(as.character(x)))

    MS1_masses = as.numeric(peaks$mz)
    MS1_RT = as.numeric(peaks$rt)/60
  }

  ####################################
  ### Extract pseudo-features with DDA
  ####################################

  if (screener == "DDA"){

    MS2_Janssen <- try(readMSData(raw_data_file, msLevel = 2,  verbose = FALSE),silent=T)
    MS1_Janssen <- try(readMSData(raw_data_file, msLevel = 1,  verbose = FALSE),silent=T)

    ### Extract useful informations from raw data:

    NS = length(MS2_Janssen) # Total number of scans
    MS2_prec_mz = precursorMz(MS2_Janssen) # Label precursor mass
    MS2_prec_rt = rtime(MS2_Janssen) # In second
    MS2_tic = sapply(1:NS,function(ttt) MS2_Janssen[[ttt]]@tic)
    polarity = polarity(MS2_Janssen)[1]

    ### Clustering

    prec_dist = proxy::dist(MS2_prec_mz, ppm_distance1)
    prec_clust = hclust(prec_dist, method = "average")
    prec_tree = cutree(prec_clust, h = ppm_search*2)
    prec_tree_class = unique(prec_tree)
    
    MS1_PEPMASS = c()
    MS1_RT = c()
    
    for (pt in prec_tree_class){
      temp_ind = which(prec_tree == pt)
      temp_mz = mean(MS2_prec_mz[temp_ind])
      temp_rt = MS2_prec_rt[temp_ind]
      
      if (length(temp_rt)>1){
        temp_rt_dist = proxy::dist(temp_rt, rt_distance1)
        temp_rt_clust = hclust(temp_rt_dist, method = "average")
        temp_rt_tree = cutree(temp_rt_clust, h = rt_search*2)
        rt_cutted = aggregate(temp_rt, list(temp_rt_tree), mean)
        MS1_PEPMASS = c(MS1_PEPMASS, rep(temp_mz, nrow(rt_cutted)))
        MS1_RT = c(MS1_RT, rt_cutted[,2]/60)
      } else {
        MS1_PEPMASS = c(MS1_PEPMASS, temp_mz)
        MS1_RT = c(MS1_RT, temp_rt/60)      }
      }

    ### Create metadata

    MS1_metadata = matrix("N/A", length(MS1_PEPMASS), 7)
    colnames(MS1_metadata) = c("PEPMASS","RT","IONMODE","ADDUCT","CHARGE","ID", "MS1_ABUNDANCE")
    MS1_metadata = data.frame(MS1_metadata)

    MS1_metadata$PEPMASS = MS1_PEPMASS
    MS1_metadata$RT = MS1_RT
    
    # Intergrate MS1 area

    mz_range = cbind(as.numeric(MS1_metadata$PEPMASS) - 0.005, as.numeric(MS1_metadata$PEPMASS) +0.005)
    rt_range = cbind(as.numeric(MS1_metadata$RT) - rt_search, as.numeric(MS1_metadata$RT) + rt_search)*60
    chrs = chromatogram(MS1_Janssen, mz = mz_range, rt = rt_range)

    for (cc in 1:NC){
      x = chrs[cc,1]@rtime
      y = chrs[cc,1]@intensity
      y[is.na(y)] = 0
      MS1_metadata$MS1_ABUNDANCE[cc] = trapz(x, y)
    }

    if (polarity==1){
      MS1_metadata$IONMODE="Positive"
      MS1_metadata$ADDUCT = "M+H"}
    if (polarity==-1){
      MS1_metadata$IONMODE="Negative"
      MS1_metadata$ADDUCT = "M-H"}

    MS1_metadata$CHARGE = 1
    MS1_metadata = MS1_metadata[order(MS1_metadata$MS1_ABUNDANCE, decreasing=T),]
    MS1_metadata = MS1_metadata[!duplicated(MS1_metadata$PEPMASS),]
    MS1_metadata$ID = paste0(basename(file_path_sans_ext(raw_data_file)),"_",1:nrow(MS1_metadata))
   }

  ##################################
  ### Filter user provided metadata
  #################################

  if (!is.null(ref)){
    ref1 = c() # New metadata
    for (i in 1:nrow(ref)){
      ppm_errors = ppm_distance1(MS1_masses,ref$PEPMASS[i])
      rt_errors =  abs(MS1_RT - as.numeric(ref$RT[i]))
      valid = which(ppm_errors<=ppm_search & rt_errors<=rt_search/60 & MS1_metadata$ADDUCT == ref$ADDUCT[i])
      nbf = length(valid)
      if (nbf>0){ # Update the found feature:
        ref_selected = do.call("rbind", replicate(nbf, ref[i,], simplify = FALSE))
        ref_selected$RT = MS1_metadata$RT[valid]
        if (nbf>1){ref_selected$ID = paste0(ref_selected$ID,"_",1:nbf)}
        ref1 = rbind.data.frame(ref1,ref_selected)
      }
    }
  } else {
    ref1 = MS1_metadata
  }
  return(ref1)
}

############################
### Internal functions:
###########################

ppm_distance1<-function(x,y){
  x = as.numeric(x)
  y = as.numeric(y)
  return(abs((x-y)/y*1000000))
}

rt_distance1<-function(x,y){
  return(abs(x-y))
}
