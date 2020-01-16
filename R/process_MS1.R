#' Read and combine targeted MS1 scans from one LC-MS/MS file
#'
#' Function used by library_generator to detect MS1 scans
#' @export

process_MS1<-function(mzdatafiles, ref, rt_search=10, ppm_search=20,
                      isomers = T,
                      baseline = 1000, relative = 5, normalized=T){

  options(stringsAsFactors = F)
  options(warn=-1)

  ### Initialize variables

  MS1_Janssen = NULL
  new_MS1_meta_data = c() # Compounds detected in the ms data
  MS1_scan_list = list() # List of spectrum2 objects
  scan_number = c() # Which scan number in the raw chromatogram
  new_PEP_mass = c() # The real mass in samples
  mass_dev = c() # Mass deviations in ppm
  spectrum_list = list() # List of spectra to save
  N=0

  ### Read the raw data file

  MS1_Janssen <- try(readMSData(mzdatafiles, msLevel = 1, verbose = FALSE),silent=T)

  if (class(MS1_Janssen)!="try-error"){ # If data contains MS1 scan

    print(paste0("Processing MS1 scans of data file ",mzdatafiles," ..."))

    ### Extract useful informations from raw data:

    NS = length(MS1_Janssen) # Total number of scans
    MS1_prec_rt = rtime(MS1_Janssen) # In second
    MS1_tic = sapply(1:NS,function(ttt) MS1_Janssen[[ttt]]@tic)
    polarity = polarity(MS1_Janssen)[1]
    if (polarity==1){ref = ref[ref$IONMODE=="Positive",]}
    if (polarity==-1){ref = ref[ref$IONMODE=="Negative",]}

    ### Extract ref:

    prec_theo=as.numeric(ref$PEPMASS)
    prec_rt=as.numeric(ref$RT)*60 # Allow N/A
    int_max_list = c() # Maximal intensity of MS1 mass

    for (i in 1:nrow(ref)){

      valid_k = 0

    # Save highest intensity MS1 data:

      if (!is.na(prec_rt[i])){
          scan_range=which(MS1_prec_rt >= prec_rt[i] - rt_search & MS1_prec_rt <= prec_rt[i] + rt_search)
     } else {scan_range=1:length(MS1_prec_rt)}

    # Find the scan that corresponds to the meta data

     if (length(scan_range)>0){
       scan_rts = c() # Validated scan retention time
       scan_tics = c()
       klist = c()

       for (k in scan_range){ # Check whether the precursor peak is detected in selected scan range
           Frag_data = MS1_Janssen[[k]]

           if (length(Frag_data@mz)>1){ # At least the scan is not empty
               ppm_dis = ppm_distance(Frag_data@mz,prec_theo[i])
               error= min(ppm_dis)
               prec_ind = which.min(ppm_dis)
               prec_int = Frag_data@intensity[prec_ind] # The intensity of precursor ion

           # We now check whether the precursor is precisely isolated
           if ((prec_int>10*baseline) & (error<ppm_search)){ # The precursor mass must be 10 times higher than baseline
               scan_rts = c(scan_rts,MS1_prec_rt[k])
               scan_tics = c(scan_tics,prec_int)
               klist = c(klist,k)
           }
      }} # End of checking scans for peaks

      if (length(klist)>0){ # Find peaks
           peaks = findpeaks(scan_tics)[,2] # intensity of precursor, find peak of XIC!
           peak_rts = scan_rts[peaks] # retention time of scans
           peak_tics = scan_tics[peaks]
           peak_range = klist[peaks]
           if (isomers){
             valid_k = separated_peaks(peak_range, peak_rts, peak_tics, rt_window = rt_search*2) # Scan number of peaks
           } else {valid_k = peak_range[which.max(peak_tics)]} # Scan number of the highest peak if no isomers
           peak_int = peak_tics[match(valid_k,peak_range)] # Precurusor intensity of selected scans
      }

    if (sum(valid_k)>0){ # If at least a scan is found

        NV = length(valid_k) # >1 if isomers present
        scan_number = c(scan_number,valid_k)  # Save scan number
        int_max_list=c(int_max_list, peak_int) # Save maximal intensity

        ### Update metadata:
        # Save detected precursor mass:

        for (vvv in valid_k){
          masslist=MS1_Janssen[[vvv]]@mz
          wp= which.min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
          dev_ppm= min(abs(masslist-prec_theo[i])/prec_theo[i]*1000000)
          dev_ppm=round(dev_ppm,2)
          new_PEP_mass = c(new_PEP_mass,masslist[wp])
          mass_dev = c(mass_dev,dev_ppm)
        }

        # Append spectra and metadata:
        for (vvv in 1:NV){
          MS1_scan_list[[N+vvv]]=MS1_Janssen[[valid_k[vvv]]]
          new_MS1_meta_data = rbind(new_MS1_meta_data,ref[i,])}
        N=N+NV
    }}} # End of the big for loop

  ### Update metadata

    new_MS1_meta_data[,"PEPMASS"]=round(as.numeric(new_PEP_mass),5)
    new_MS1_meta_data[,"RT"]= round(MS1_prec_rt[scan_number]/60,2)
    new_MS1_meta_data[,"FILENAME"]=rep(basename(mzdatafiles),N)
    new_MS1_meta_data[,"MSLEVEL"]=rep(1,N)
    new_MS1_meta_data[,"TIC"]= int_max_list
    new_MS1_meta_data[,"PEPMASS_DEV"]=mass_dev
    new_MS1_meta_data[,"SCAN_NUMBER"] = scan_number

  ### Update metadata with library search parameters

    new_MS1_meta_data[,"PARAM_RT_SEARCH"]= rep(rt_search,N)
    new_MS1_meta_data[,"PARAM_MASS_SEARCH_PPM"]= rep(ppm_search,N)
    new_MS1_meta_data[,"PARAM_BASELINE_INTENSITY"]= rep(baseline,N)
    new_MS1_meta_data[,"PARAM_RELATIVE_INTENSITY"]= rep(relative,N)
    if (normalized){new_MS1_meta_data[,"PARAM_NORMALIZED"]= rep("Yes",N)
    } else {new_MS1_meta_data[,"PARAM_NORMALIZED"]= rep("No",N)}

  ### Denoise spectra

  if (!is.null(new_MS1_meta_data)){
    included=c()
    n0=0

    for (i in 1:N){
      dat = cbind(MS1_scan_list[[i]]@mz,MS1_scan_list[[i]]@intensity)
      # Filter background noise and choose only peaks until +7Da of the precursor!!

      baseline1= max(baseline,max(dat[,2])*relative/100)
      selected = which(dat[,1]>new_MS1_meta_data$PEPMASS[i]-0.5 & dat[,1] < new_MS1_meta_data$PEPMASS[i]+7 & dat[,2]>baseline1)

      if (length(selected)>1){ # At least 2 peaks should be present for isotope determination
        dat = dat[selected,]
        dat = matrix(dat,ncol=2)
        if (normalized){dat[,2]=dat[,2]/max(dat[,2])*100}
        n0=n0+1
        spectrum_list[[n0]]=dat
        included= c(included,i)}}

    ### Keep only no empty spectra
    new_MS1_meta_data = new_MS1_meta_data[included,]
  }
  } else {
   print(paste0("No MS1 scan in the data file ",mzdatafiles," !"))
 }

  return(list(sp=spectrum_list,metadata=new_MS1_meta_data))
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

