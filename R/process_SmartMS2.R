#' Read and combine targeted MS2 scans from one LC-MS/MS file using default MergeION algorithm
#'
#' Function used by library_generator to detect MS2 scans
#' 
#' @importFrom BiocGenerics basename
#' @importFrom MSnbase readMSData rtime tic fData readMgfData precursorMz polarity chromatogram
#' @importFrom pracma findpeaks trapz
#' @importFrom mzR openMSfile header peaks
#'
#' @export
#' 
process_SmartMS2<-function(mzdatafiles = NULL, ref = NULL, 
      rt_search = 10, rt_gap = 30, ppm_search = 10, mz_search = 0.01,
      baseline= 1000, relative = 5, max_peaks = 200, normalized=T){

 options(stringsAsFactors = F)
 options(warn=-1)
 
 if ("FILENAME" %in% colnames(ref)){
   valid = c(which(basename(ref$FILENAME) == basename(mzdatafiles)),  which(ref$FILENAME =="N/A"))
   ref = ref[valid,,drop=FALSE]
 }

 ### Initialize variables

 new_MS2_meta_data = c() # Compounds detected in the ms data
 MS2_scan_list = list() # List of spectrum2 objects
 scan_number = c() # Which scan number in the raw chromatogram is found
 new_PEP_mass = c() # The real mass in samples
 mass_dev = c() # Mass deviations in ppm
 spectrum_list = list() # List of spectra to save
 NNN=0 # Numerator

 #####################
 ### Load MS2 Scans###
 #####################
 if (nrow(ref)==0){MS2_Janssen=NULL}else{
 
  MS2_Janssen <- try(readMSData(mzdatafiles, msLevel = 2, verbose = FALSE, mode = "inMemory",  centroided = T),silent=T)
 
 if (class(MS2_Janssen)=="try-error"){MS2_Janssen=NULL}
 }

 
 if (!is.null(MS2_Janssen)){ # If data contains MS2 scan

   print("Processing MS2 scans of data file ...")
   
   ### Extract useful informations from raw data:

   NS = length(MS2_Janssen) # Total number of scans
   MS2_prec_mz = precursorMz(MS2_Janssen) # Label precursor mass
   targets =  unique(MS2_prec_mz) # All targeted masses
   MS2_prec_rt = rtime(MS2_Janssen) # In second
   MS2_tic = as.numeric(MSnbase::tic(MS2_Janssen))

   ### Filter ref because not all targeted m/z exists or fragmented in the sample! Important for not to search whatever

   dev_targets = sapply(as.numeric(ref$PEPMASS),function(x) min(abs(x-targets)))
   valid = which(dev_targets <= 1) #  Find targeted metadata in experimental file!! - faster!

   ref = ref[valid,,drop=FALSE]
   prec_theo= as.numeric(ref$PEPMASS)
   prec_rt=as.numeric(ref$RT)*60 # Allow N/A
   
   ####################################################
   ### Go through each ref item to find adequate scan##
   ####################################################
   
   if (nrow(ref)>0){

      for (i in 1:nrow(ref)){
        
        # 1. Define a wide search range based on targeted precursor mass
        
        scan_range = which(abs(MS2_prec_mz-prec_theo[i])<1 & MS2_tic>baseline) 
        if (!is.na(prec_rt[i])){
            time_range = which(MS2_prec_rt >= prec_rt[i] - rt_search & MS2_prec_rt <= prec_rt[i] + rt_search)
            scan_range = intersect(scan_range,time_range)
        }
        
        scan_range = sort(unique(as.numeric(scan_range)))
        
        # 2. Refine the search range based on exact precursor mass or mass detected in actual spectrum
        
        if (length(scan_range)>0){
          
          temp_scan_range = c() # Validated scan number
          temp_mz0 = c()
          
          for (k in scan_range){
            
            checked = 0 # The scan validation state
            Frag_data = MS2_Janssen[[k]]
            
            ppm_error = ppm_distance(MS2_prec_mz[k], prec_theo[i])
            abs_error = abs(MS2_prec_mz[k]-prec_theo[i])
            
            # Easy situation if scans are labeled correctly:
            
            if (ppm_error<=ppm_search || abs_error<=mz_search){
              checked = 1
              mz0 = MS2_prec_mz[k]
            }
          
            # If the scan is not correctly labeled, check in details actual spectrum:
            
            if (checked==0  & length(Frag_data@mz)>1){
              ppm_dis = ppm_distance(Frag_data@mz,prec_theo[i])
              mz_dev = abs(Frag_data@mz - prec_theo[i])
              prec_ind = which.min(ppm_dis)
              prec_int = Frag_data@intensity[prec_ind] # The intensity of precursor ion

              if (prec_int>baseline*2 & (ppm_dis[prec_ind]<=ppm_search || mz_dev[prec_ind]<=mz_search)){
                checked = 1
                mz0 = Frag_data@mz[prec_ind]
              } # The precursor mass must be higher than baseline
            }
            
            # Add to filters scans  
              
            if (checked==1){
              temp_scan_range = c(temp_scan_range, k)
              temp_mz0 = c(temp_mz0, mz0)
            }
          }
          
          scan_range = temp_scan_range
          scan_mz0 = temp_mz0
          
      ### 3. Select "best" scan and calculate deviation
    
          if (length(scan_range)>0){
            
            scan_rts = MS2_prec_rt[scan_range]
            scan_tics = MS2_tic[scan_range]
            ind_features = separated_peaks2(scan_range, scan_rts, scan_tics, rt_gap)

            valid_k = sapply(ind_features, function(x) x[1]) # Highest TIC scan

            NV = length(valid_k) # >1 if isomers present
            mz = scan_mz0[match(valid_k,scan_range)] # Precursor m/z of valid scan
            dev_ppm= round(sapply(mz, function(x) min(ppm_distance(x,prec_theo[i]))),2)

            scan_number = c(scan_number,valid_k)  # Save scan numbers
            new_PEP_mass = c(new_PEP_mass, mz)
            mass_dev = c(mass_dev,dev_ppm)
        
          ### 4. Calculate smart ion append spectra and metadata:
      
            for (vvv in 1:NV){
              
                ind_feature = as.numeric(ind_features[[vvv]])
                tmp_scans = MS2_Janssen[ind_feature]

                scan_best_tic = MS2_Janssen[[ind_feature[as.numeric(which.max(MSnbase::tic(tmp_scans)))]]]
                scan_best_tic = cbind(scan_best_tic@mz, scan_best_tic@intensity)

                scan_most_peaks = MS2_Janssen[[ind_feature[as.numeric(which.max(mzR::peaksCount(tmp_scans)))]]]
                scan_most_peaks = cbind(scan_most_peaks@mz, scan_most_peaks@intensity)

                scan_least_peaks = MS2_Janssen[[ind_feature[as.numeric(which.min(mzR::peaksCount(tmp_scans)))]]]
                scan_least_peaks = cbind(scan_least_peaks@mz, scan_least_peaks@intensity)

                scan_all = list(scan_best_tic, scan_most_peaks, scan_least_peaks)
                
                scan_merged = average_spectrum(scan_all, mz_window = mz_search*2)
                scan_merged = scan_merged$new_spectrum
                
                MS2_scan_list[[NNN+vvv]] = scan_merged
                new_MS2_meta_data = rbind(new_MS2_meta_data,ref[i,,drop=FALSE])
            }
            
       NNN=NNN+NV
      }
    }
  } # End of screening all reference masses
  
  #####################   
  ### Create metadata##
  #####################
    
  if (NNN>0){ 
    
    new_MS2_meta_data[,"PEPMASS"] = round(as.numeric(new_PEP_mass),5)
    new_MS2_meta_data[,"RT"] = round(MS2_prec_rt[scan_number]/60,2)  # minutes
    new_MS2_meta_data[,"FILENAME"] = rep(basename(mzdatafiles),NNN)
    new_MS2_meta_data[,"MSLEVEL"] = rep(2,NNN)
    new_MS2_meta_data[,"TIC"] = MS2_tic[scan_number]
    new_MS2_meta_data[,"PEPMASS_DEV"] = mass_dev
    new_MS2_meta_data[,"SCAN_NUMBER"] = scan_number
    
    #################################
    ### Process and collect spectra##
    #################################
  
    included = c() # Not filtered
    n0=0
      
    for (i in 1:NNN){
        
       sp0 = MS2_scan_list[[i]]
       
       if (!is.null(sp0)){
         
          if (nrow(sp0)>2){
            sp1 = denoise_ms2_spectrum(sp0, new_MS2_meta_data$PEPMASS[i], max_peaks, relative, normalized)
        
            if (!is.null(sp1)){
          
            if (nrow(sp1)>1){
              included = c(included, i)
              n0 = n0 + 1
              spectrum_list[[n0]]=sp1
        }}}}
    }
    
    new_MS2_meta_data = new_MS2_meta_data[included,,drop=FALSE]
    id_kept = unique(new_MS2_meta_data$ID)
    ref = ref[match(id_kept,ref$ID),,drop=FALSE]
   } 
   }}

  if (!is.null(new_MS2_meta_data)){
    if (nrow(new_MS2_meta_data)==0){  
      print(paste0("No MS2 scan in the data file ",mzdatafiles," matches with metadata!"))
  }} else { print(paste0("No MS2 scan in the data file ",mzdatafiles," matches with metadata!"))}

  return(list(sp=spectrum_list,metadata=new_MS2_meta_data,ref_MS2=ref))
}

###########################
### Internal functions:####
###########################

#' ppm error calculation:
#' @param x
#' @param y

ppm_distance<-function(x,y){
  x = as.numeric(x)
  y = as.numeric(y)
  if (y>100){
    ppm = abs((x-y))/y*1000000
  } else {
    ppm = abs(x-y)
    ppm[ppm<0.01]=0
  }
  return(ppm)
}

#' Find indexes from "ranges": separated isomer peaks according to rt and intensity
#'@param ranges
#'@param rts
#'@param tics
#'@param rt_gap

separated_peaks2<-function(ranges, rts, tics, rt_gap){

  # ranges: scan or peak number
  NR = length(ranges)
  ind_features = list()
  kf = 0
  
 # screen:
  if (NR>1){
    tmp = data.matrix(cbind(ranges,rts,tics))
    tmp = tmp[order(tmp[,3],decreasing=T),] # Order the tmp by intensity
    while (nrow(tmp)>0){
      kf = kf +1
      feature_rt = tmp[1,2]
      valid = which(abs(tmp[,2]-feature_rt)<=rt_gap)
      ind_features[[kf]] = tmp[valid, 1]
      tmp  = tmp[-valid,,drop = FALSE]
  }} else {ind_features[[1]] = ranges}
  
  # check:
  return(ind_features)
}

#' Keep top peaks
#'@param sp
#'@param mz0
#'@param max_peak
#'@param min_relative
#'@param normalized

denoise_ms2_spectrum<-function(sp, mz0, max_peak, min_relative, normalized = T){
  
  denoised_spectrum = matrix(c(0,0),1,2)
  
  if (nrow(sp)>0){
    
    # Check resolution:
    
    checked = any(sapply(sp[,1], decimalplaces)>2) # At least 2 values after decimal
    
    # Filter top peaks:
    
    sp = sp[order(sp[,2], decreasing = T),,drop=FALSE]
    tops = min(max_peak, nrow(sp))  
    sp = sp[1:tops,,drop=FALSE]
    
    # Normalize to 100:
    
    sp1 = sp
    sp1[,2] = sp1[,2]/max(sp1[,2])*100
    
    # Relative Intensity filter:
    
    #filter = which(sp1[,2]>=min_relative & sp1[,1]<mz0-1)
    ## keep the PEPMASS in the final spectra libary 
    filter = which(sp1[,2]>=min_relative & sp1[,1] <= mz0)
    if (normalized){sp = sp1}  
    sp = sp[filter,,drop=FALSE]
    
    # Check validity:
    
    if (nrow(sp)>0 & checked){
      sp = sp[order(sp[,1]),,drop=FALSE]
      if (normalized){sp[,2] = sp[,2]/max(sp[,2])*100}
      denoised_spectrum = sp
    }
  }
  return(denoised_spectrum)
}

decimalplaces <- function(x){
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
