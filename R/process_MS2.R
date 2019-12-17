#' Read and combine targeted MS2 scans from one LC-MS/MS file
#'
#' Function used by library_generator to detect MS2 scans
#' @export

process_MS2<-function(mzdatafiles, ref, rt_search=10, ppm_search=20,
                      MS2_type = c("DDA","Targeted"), isomers = T,
                      baseline= 1000, relative = 5,normalized=T){

 options(stringsAsFactors = F)
 options(warn=-1)

 ### Initialize variables

 MS2_Janssen = NULL
 new_MS2_meta_data = c() # Compounds detected in the ms data
 MS2_scan_list = list() # List of spectrum2 objects
 scan_number = c() # Which scan number in the raw chromatogram is found
 new_PEP_mass = c() # The real mass in samples
 mass_dev = c() # Mass deviations in ppm
 spectrum_list = list() # List of spectra to save
 N=0 # Numerator

### Read the raw data file

 MS2_Janssen <- try(readMSData(mzdatafiles, msLevel = 2,  verbose = FALSE),silent=T)

 if (class(MS2_Janssen)!="try-error"){ # If data contains MS2 scan

    print(paste0("Processing MS2 scans of data file ",mzdatafiles," ..."))

   ### Extract useful informations from raw data:

    NS = length(MS2_Janssen) # Total number of scans
    MS2_prec_mz = precursorMz(MS2_Janssen) # Label precursor mass
    targets =  unique(MS2_prec_mz) # All targeted masses
    MS2_prec_rt = rtime(MS2_Janssen) # In second
    MS2_tic = sapply(1:NS,function(ttt) MS2_Janssen[[ttt]]@tic)
    polarity = polarity(MS2_Janssen)[1]

  ### Filter ref because not all targeted m/z exists or fragmented in the sample! Important for not to search whatever

    dev_targets = sapply(ref$PEPMASS,function(x) min(abs(x-targets)))
    valid = which(dev_targets <= 1) #  Find targeted metadata in experimental file!! - faster!
    if (polarity==1){valid = intersect(valid, which(ref$IONMODE=="Positive"))}
    if (polarity==-1){valid = intersect(valid, which(ref$IONMODE=="Negative"))}

    ref = ref[valid,]
    prec_theo=ref$PEPMASS
    prec_rt=as.numeric(ref$RT)*60 # Allow N/A

    if (nrow(ref)>0){

    ### Check one by one targeted m/z:

      for (i in 1:nrow(ref)){

        valid_k = 0

     # Define search range - precursor mass + RT if precised:

        if (MS2_type=="Targeted"){scan_range = which(abs(MS2_prec_mz-prec_theo[i])<1)}
        if (MS2_type=="DDA"){scan_range = which(ppm_distance(MS2_prec_mz,prec_theo[i])<ppm_search)}
        scan_range = as.numeric(scan_range)
        if (!is.na(prec_rt[i])){
          time_range = which(MS2_prec_rt >= prec_rt[i] - rt_search & MS2_prec_rt <= prec_rt[i] + rt_search)
          scan_range = intersect(scan_range,time_range)
        }

       if (length(scan_range)<2){}
       if (length(scan_range)>=1){

    ### Targeted mode: check inside scans the precursor masses:

        if (MS2_type=="Targeted"){
          scan_rts = c() # Validated scan retention time
          scan_tics = c()
          klist = c()
          for (k in scan_range){
            Frag_data = MS2_Janssen[[k]]
            if (length(Frag_data@mz)>1){
               ppm_dis = ppm_distance(Frag_data@mz,prec_theo[i])
               error= min(ppm_dis)
               prec_ind = which.min(ppm_dis)
               prec_int = Frag_data@intensity[prec_ind] # The intensity of precursor ion
               if ((prec_int>baseline) & (error<ppm_search)){ # The precursor mass must be higher than baseline
                 scan_rts = c(scan_rts,MS2_prec_rt[k])
                 scan_tics = c(scan_tics,MS2_tic[k])
                 klist = c(klist,k)
            }}
          } # End of checking scans for peaks
        if (length(klist)>0){ # Find peaks
            peaks = findpeaks(scan_tics)[,2]
            peak_rts = scan_rts[peaks] # retention time of scans
            peak_tics = scan_tics[peaks]
            peak_range = klist[peaks]
            if (isomers & is.na(prec_rt[i])){ # Report several peaks only if retention time is not precised
                valid_k = separated_peaks(peak_range, peak_rts, peak_tics, rt_window = rt_search*2)
          } else {valid_k = peak_range[which.max(peak_tics)]}
      }}

    ### DDA mode: simply check scan labels in the scan range:

      if (MS2_type=="DDA"){
        scan_tics = MS2_tic[scan_range] # peaks = targeted scans
        scan_rts = MS2_prec_rt[scan_range]
        if (isomers){
        valid_k = separated_peaks(scan_range, scan_rts, scan_tics, rt_window = rt_search*2)
        } else {valid_k = scan_range[which.max(scan_tics)]}
      }

    # When search finished:

    if (sum(valid_k)==0){} # If no scan ever found do nothing

    if (sum(valid_k)!=0){ # If at least one scan is found

      NV = length(valid_k) # >1 if isomers present
      scan_number = c(scan_number,valid_k)  # Save scan number

      ### Update metadata:
      # Save detected precursor mass:

      if (MS2_type=="DDA"){
          mz = MS2_prec_mz[valid_k]
          dev_ppm= sapply(mz, function(x) min(ppm_distance(x,prec_theo[i])))
          dev_ppm=round(dev_ppm,2)
          new_PEP_mass = c(new_PEP_mass,mz)
          mass_dev = c(mass_dev,dev_ppm)
        }

      if (MS2_type=="Targeted"){
          for (vvv in valid_k){
            masslist=MS2_Janssen[[vvv]]@mz
            wp= which.min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
            dev_ppm= min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
            dev_ppm=round(dev_ppm,2)
            new_PEP_mass = c(new_PEP_mass,masslist[wp])
            mass_dev = c(mass_dev,dev_ppm)
          }
        }

    # Append spectra and metadata:
      for (vvv in 1:NV){
        MS2_scan_list[[N+vvv]]=MS2_Janssen[[valid_k[vvv]]]
        new_MS2_meta_data = rbind(new_MS2_meta_data,ref[i,])}
      N=N+NV
    }}} # End of big for loop of scan matching to reference masses

  ### Update metadata
    new_MS2_meta_data[,"PEPMASS"] = round(as.numeric(new_PEP_mass),5)
    new_MS2_meta_data[,"RT"] = round(MS2_prec_rt[scan_number]/60,2)  # minutes
    new_MS2_meta_data[,"FILENAME"] = rep(basename(mzdatafiles),N)
    new_MS2_meta_data[,"MSLEVEL"] = rep(2,N)
    new_MS2_meta_data[,"TIC"] = MS2_tic[scan_number]
    new_MS2_meta_data[,"PEPMASS_DEV"] = mass_dev
    new_MS2_meta_data[,"SCAN_NUMBER"] = scan_number

    ### Update metadata with library search parameters

    new_MS2_meta_data[,"PARAM_RT_SEARCH"]= rep(rt_search,N)
    new_MS2_meta_data[,"PARAM_MASS_SEARCH_PPM"]= rep(ppm_search,N)
    new_MS2_meta_data[,"PARAM_BASELINE_INTENSITY"]= rep(baseline,N)
    new_MS2_meta_data[,"PARAM_RELATIVE_INTENSITY"]= rep(relative,N)
    if (normalized){new_MS2_meta_data[,"PARAM_NORMALIZED"]= rep("Yes",N)
    } else {new_MS2_meta_data[,"PARAM_NORMALIZED"]= rep("No",N)}

  ### Denoise spectra

  if (!is.null(new_MS2_meta_data)){
    included = c() # Not filtered
    n0=0

    for (i in 1:N){

      dat = cbind(MS2_scan_list[[i]]@mz,MS2_scan_list[[i]]@intensity)
      baseline1= max(baseline,max(dat[,2])*relative/100)

      # Cut only masses smaller than precursor and filter background noise:
      selected = which((dat[,1] < new_MS2_meta_data$PEPMASS[i]+10) & (dat[,2]>baseline1))

      if (length(selected)>0){
        dat = dat[selected,]
        dat = matrix(dat,ncol=2)
        if (normalized){dat[,2]=dat[,2]/max(dat[,2])*100}
        n0=n0+1
        spectrum_list[[n0]]=dat
        included= c(included,i)}}

  ### Keep only no empty spectra
      new_MS2_meta_data = new_MS2_meta_data[included,]
      id_kept = unique(new_MS2_meta_data$ID)
      ref = ref[match(id_kept,ref$ID),]
     }
    } else {
      print(paste0("No MS2 scan in the data file ",mzdatafiles," matches with metadata!"))}
   } else {
  print(paste0("No MS2 scan in the data file ",mzdatafiles," !"))
 }

  return(list(sp=spectrum_list,metadata=new_MS2_meta_data,ref_MS2=ref))
}

############################
### Internal functions:
###########################

ppm_distance<-function(x,y){
  return(abs((x-y)/y*1000000))}

# Find indexes of separated peaks according to rt and intensity
separated_peaks<-function(ranges, rts, tics, rt_window=20){

  # ranges: scan or peak number
  NR = length(ranges)

 # screen:
  if (NR>1){
    tmp = data.matrix(cbind(ranges,rts,tics))
    tmp = tmp[order(tics,decreasing=T),] # Order the tmp by intensity
    valid_k = tmp[1,1]
    previous_rt = tmp[1,2]
    for (i in 2:NR){
       dis = abs(tmp[i,2]-previous_rt)
       valid = which(dis<=rt_window) # Check peak overlap
       if (length(valid)==0){ # No overlap
         valid_k = c(valid_k,tmp[i,1])
         previous_rt = c(previous_rt,tmp[i,2])}
   }} else {valid_k = ranges}
  return(unique(valid_k))
}

