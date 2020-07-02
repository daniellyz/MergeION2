#' Read and combine targeted MS1 and MS2 scans from one LC-MS/MS file using RMassBank algorithm
#'
#' Function used by library_generator
#' 
#' @importFrom BiocGenerics basename
#' @importFrom mzR openMSfile close peaks
#' 
#' @export
#' 
process_RMassBank<-function(mzdatafiles = NULL, ref = NULL, polarity = c("Positive", "Negative"), include.MS1 = FALSE,
                            add.adduct = FALSE, rt_search = 10, ppm_search = 10, 
                            baseline= 1000, relative = 5, max_peaks = 200, recalibration = F, normalized=T){
  
  options(stringsAsFactors = F)
  options(warn=-1)
  library(RMassBank)
  
  # @importFrom RMassBank ppm RmbSettingsTemplate loadRmbSettings loadList newMsmsWorkspace analyzeMsMs aggregateSpectra makeRecalibration recalibrateSpectra findLevel findMz findRt findName findMsMsHRperxcms findFormula archiveResults cleanElnoise reanalyzeFailpeaks filterMultiplicity processProblematicPeaks deprofile.scan
  
  if (add.adduct){mode = c("pH", "pNa", "pM")}
  if (!add.adduct){mode = c("pH")}
  
  if (polarity=="Positive"){ref = ref[ref$ADDUCT=="M+H",,drop=FALSE]}
  if (polarity=="Negative"){ref = ref[ref$ADDUCT=="M-H",,drop=FALSE]}
  
  output_metadata = c()
  output_sp = list()
  kkk = 0

  ###########################
  ######Prepare metadata#####
  ###########################
  
  valid = which(basename(ref$FILENAME) == basename(mzdatafiles)) 

  if (length(valid)>0){
  
    ref = ref[valid,,drop=FALSE]
    ref0 = ref # Backup!
    
    ref$ID = 1:nrow(ref)
    
    if (!("Name" %in% colnames(ref))){ref$Name = ref$ID}
    if (!("CAS" %in% colnames(ref))){ref$CAS = ref$ID}
  
    #############################
    ######Prepare setting file####
    ##############################
  
    RmbSettingsTemplate("mysettings.ini")
    template = readLines("mysettings.ini")

    ind_rt_search = grep('rtMargin:',template)
    template[ind_rt_search] = paste0("rtMargin: ", rt_search/60)
  
    ind_rt_search = grep('rtShift:',template)
    template[ind_rt_search] = paste0("rtShift: ", 0)
  
    ind_ppm_search = grep('ppmHighMass:',template)
    template[ind_ppm_search] = paste0("    ppmHighMass: ", ppm_search)
  
    ind_ppm_search2 = grep('ppmFine:',template)[1]
    template[ind_ppm_search2] = paste0("    ppmFine: ", ppm_search/2)
  
    ind_baseline = grep('prelimCut:',template)
    template[ind_baseline] = paste0("    prelimCut: ", baseline)
  
    ind_ppm_search3 = grep('ppmFine:',template)[2]
    template[ind_ppm_search3] = paste0("    ppmFine: ", ppm_search)
  
    writeLines(template, con = "mysettings.ini")
  
    #####################
    ######Processing#####
    #####################
  
    for (rr in 1:nrow(ref)){
      
      ref2 = ref[rr,,drop=FALSE]
      for (x in 1:ncol(ref2)){ref2[1,x] = gsub(",","_", ref2[1,x])} # Pb with compound list...
      
      ref2$PEPMASS = as.numeric(ref2$PEPMASS)
      write.csv(ref2, "Compoundlist.csv", sep=",", col.names = T, row.names = F, quote = F)
      loadRmbSettings("mysettings.ini")
      loadList("Compoundlist.csv")
      
      w <- newMsmsWorkspace()
      w@files = mzdatafiles
      
      # Basic processing:
      w1 <- msmsWorkflow1(w = w, cpdid = ref2$ID, mode = mode, steps = c(1:2), readMethod = "mzR")
      sp1 = w1@spectra[[1]]@parent
      sp1 = cbind(sp1@mz, sp1@intensity)
      temp_sp2 = w1@spectra[[1]]@children # MS2
      
      prec_mz = w1@spectra[[1]]@mz
      prec_rt_ms1 = round(w1@spectra[[1]]@parent@rt/60,2)
      prec_rt_ms2 = round(median(sapply(temp_sp2, function(x) x@rt))/60, 2)
      prec_scan_ms1 = w1@spectra[[1]]@parent@acquisitionNum
      prec_scan_ms2 =  round(median(sapply(temp_sp2, function(x) x@acquisitionNum)),0)
      
      # Aggregation and recalibration:
      
      if (length(temp_sp2)>1){ # Further aggregate if multiple scans found
        w1 <- msmsWorkflow1(w = w1, cpdid = ref2$ID, mode = mode, steps = 3, readMethod = "mzR")
        sp2 = cbind(w1@aggregated$mzFound, w1@aggregated$intensity)
        if (recalibration){
            w1 <- msmsWorkflow1(w = w1, cpdid = ref2$ID, mode = mode, steps = c(4:6), readMethod = "mzR")
            sp2 = cbind(w1@aggregated$mzFound, w1@aggregated$intensity)
            sp1 = cbind(w1@spectra[[1]]@parent@mz, w1@spectra[[1]]@parent@intensity)  
          }
      } 
      
      if (length(temp_sp2)==1){sp2 = cbind(temp_sp2[[1]]@mz, temp_sp2[[1]]@intensity)}
      
      # Add MS2     

      if (w1@spectra[[1]]@found){
        
        temp_metadata_ms2 = ref2
        temp_metadata_ms2$PEPMASS = prec_mz 
        temp_metadata_ms2$RT = prec_rt_ms2
        temp_metadata_ms2$MSLEVEL=2  
        temp_metadata_ms2$TIC = sum(sp2[,2])
      
        temp_metadata_ms2$PEPMASS_DEV = (prec_mz - ref2$PEPMASS)/ref2$PEPMASS*1000000
        temp_metadata_ms2$SCAN_NUMBER = prec_scan_ms2 # Median scan of aggregated
        temp_metadata_ms2$ID = ref0$ID[rr]
        
        temp_sp_ms2 = denoise_ms2_spectrum(sp2, prec_mz, max_peaks, relative, normalized)
        
        kkk = kkk+1
        output_metadata = rbind(output_metadata, temp_metadata_ms2)
        output_sp[[kkk]] = temp_sp_ms2
        
        # Add MS1
        
        if (include.MS1){
          
          temp_metadata_ms1 = ref2
          temp_metadata_ms1$PEPMASS = prec_mz 
          temp_metadata_ms1$RT = prec_rt_ms1
          temp_metadata_ms1$MSLEVEL= 1
          temp_metadata_ms1$TIC = sum(sp1[,2])
          
          temp_metadata_ms1$PEPMASS_DEV = (prec_mz - ref2$PEPMASS)/ref2$PEPMASS*1000000
          temp_metadata_ms1$SCAN_NUMBER = prec_scan_ms1 # Median scan of aggregated
          temp_metadata_ms1$ID = ref0$ID[rr]
          
          temp_sp_ms1 = denoise_ms1_spectrum(sp1, prec_mz, max_peaks, relative, normalized)
          
          kkk = kkk+1
          output_metadata = rbind(output_metadata, temp_metadata_ms1)
          output_sp[[kkk]] = temp_sp_ms1
        }
      }
    }

   ###########################
   ####Filter empty spectra###
   ###########################
      
  filter = which(sapply(output_sp, nrow)>0)
  output_metadata = output_metadata[filter,,drop=FALSE]
  output_sp = output_sp[filter]
  }
  
  if (!is.null(dev.list())){dev.off()}
  return(list(sp=output_sp,metadata=output_metadata))
}


#############################
######Internal Functions#####
#############################

msmsWorkflow1<-function(w, cpdid, mode = "pH", steps = c(1:8), confirmMode = FALSE, 
                         newRecalibration = TRUE, useRtLimit = TRUE, archivename = NA, 
                         readMethod = "mzR", findPeaksArgs = NULL, plots = FALSE, 
                         precursorscan.cf = FALSE, settings = getOption("RMassBank"), 
                         analyzeMethod = "formula", progressbar = "progressBarHook", 
                         MSe = FALSE){
  
  if (!any(mode %in% c("pH", "pNa", "pNH4", 
                       "pM", "mH", "mFA", "mM", ""))) 
    stop(paste("The ionization mode", mode, "is unknown."))
  if (!is.na(archivename)) 
    w@archivename <- archivename
  nProg <- 0
  nLen <- length(w@files)
  allUnknown <- FALSE
  #if (all(.listEnvEnv$listEnv$compoundList$Level == "5")) {
  #  allUnknown <- TRUE
  #  message("All compounds are unknown, the workflow will be adjusted accordingly")
  #}
  if (readMethod == "minimal") {
    opt <- getOption("RMassBank")
    opt$recalibrator$MS1 <- "recalibrate.identity"
    opt$recalibrator$MS2 <- "recalibrate.identity"
    opt$add_annotation <- FALSE
    opt$multiplicityFilter <- 1
    options(RMassBank = opt)
    settings <- getOption("RMassBank")
    analyzeMethod <- "intensity"
  }
  if (!all(steps > 4) & !is.null(w@parent)) {
    rc <- w@rc
    rc.ms1 <- w@rc.ms1
    w <- w@parent
    w@rc <- rc
    w@rc.ms1 <- rc.ms1
  }
  if (1 %in% steps) {
    message("msmsWorkflow: Step 1. Acquire all MSMS spectra from files")
    w <- msmsRead1(w = w, files = w@files, cpdid = cpdid, readMethod = readMethod, 
                  mode = mode, confirmMode = confirmMode, useRtLimit = useRtLimit, 
                  Args = findPeaksArgs, settings = settings, progressbar = progressbar, 
                  MSe = MSe)
  }
  if (2 %in% steps) {
    nProg <- 0
    message("msmsWorkflow: Step 2. First analysis pre recalibration")
    if (allUnknown) {
      analyzeMethod <- "intensity"
    }
    pb <- do.call(progressbar, list(object = NULL, value = 0, 
                                    min = 0, max = nLen))
    w@spectra <- as(lapply(w@spectra, function(spec) {
      s <- analyzeMsMs(spec, mode = mode, detail = TRUE, 
                       run = "preliminary", filterSettings = settings$filterSettings, 
                       spectraList = settings$spectraList, method = analyzeMethod)
      nProg <<- nProg + 1
      pb <- do.call(progressbar, list(object = pb, value = nProg))
      return(s)
    }), "SimpleList")
    suppressWarnings(do.call(progressbar, list(object = pb, 
                                               close = TRUE)))
  }
  if (3 %in% steps) {
    message("msmsWorkflow: Step 3. Aggregate all spectra")
    w@aggregated <- aggregateSpectra(w@spectra, addIncomplete = TRUE)
  }
  if (allUnknown) {
    w@aggregated$noise <- FALSE
    w@aggregated$noise <- FALSE
    w@aggregated$reanalyzed.formula <- NA
    w@aggregated$reanalyzed.mzCalc <- NA
    w@aggregated$reanalyzed.dppm <- NA
    w@aggregated$reanalyzed.formulaCount <- NA
    w@aggregated$reanalyzed.dbe <- NA
    w@aggregated$matchedReanalysis <- NA
    w@aggregated$filterOK <- TRUE
    w@aggregated$problematicPeak <- FALSE
    w@aggregated$formulaMultiplicity <- unlist(sapply(table(w@aggregated$cpdID), 
                                                      function(x) rep(x, x)))
    return(w)
  }
  if (4 %in% steps) {
    message("msmsWorkflow: Step 4. Recalibrate m/z values in raw spectra")
    if (newRecalibration) {
      recal <- makeRecalibration(w, mode, recalibrateBy = settings$recalibrateBy, 
                                 recalibrateMS1 = settings$recalibrateMS1, recalibrator = settings$recalibrator, 
                                 recalibrateMS1Window = settings$recalibrateMS1Window)
      w@rc <- recal$rc
      w@rc.ms1 <- recal$rc.ms1
    }
    w@parent <- w
    w@aggregated <- data.frame()
    spectra <- recalibrateSpectra(mode, w@spectra, w = w, 
                                  recalibrateBy = settings$recalibrateBy, recalibrateMS1 = settings$recalibrateMS1)
    w@spectra <- spectra
  }
  if (5 %in% steps) {
    nProg <- 0
    message("msmsWorkflow: Step 5. Reanalyze recalibrated spectra")
    pb <- do.call(progressbar, list(object = NULL, value = 0, 
                                    min = 0, max = nLen))
    w@spectra <- as(lapply(w@spectra, function(spec) {
      if (findLevel(spec@id, TRUE) == "unknown") {
        analyzeMethod <- "intensity"
      }
      else {
        analyzeMethod <- "formula"
      }
      s <- analyzeMsMs(spec, mode = mode, detail = TRUE, 
                       run = "recalibrated", filterSettings = settings$filterSettings, 
                       spectraList = settings$spectraList, method = analyzeMethod)
      nProg <<- nProg + 1
      pb <- do.call(progressbar, list(object = pb, value = nProg))
      return(s)
    }), "SimpleList")
    suppressWarnings(do.call(progressbar, list(object = pb, 
                                               close = TRUE)))
    do.call(progressbar, list(object = pb, close = TRUE))
  }
  if (6 %in% steps) {
    message("msmsWorkflow: Step 6. Aggregate recalibrated results")
    w@aggregated <- aggregateSpectra(w@spectra, addIncomplete = TRUE)
    if (!is.na(archivename)) 
      archiveResults(w, paste(archivename, ".RData", 
                              sep = ""), settings)
    w@aggregated <- cleanElnoise(w@aggregated, settings$electronicNoise, 
                                 settings$electronicNoiseWidth)
  }
  if (7 %in% steps) {
    message("msmsWorkflow: Step 7. Reanalyze fail peaks for N2 + O")
    w@aggregated <- reanalyzeFailpeaks(w@aggregated, custom_additions = "N2O", 
                                       mode = mode, filterSettings = settings$filterSettings, 
                                       progressbar = progressbar)
    if (!is.na(archivename)) 
      archiveResults(w, paste(archivename, "_RA.RData", 
                              sep = ""), settings)
  }
  if (8 %in% steps) {
    message("msmsWorkflow: Step 8. Peak multiplicity filtering")
    if (is.null(settings$multiplicityFilter)) {
      message("msmsWorkflow: Step 8. Peak multiplicity filtering skipped because multiplicityFilter parameter is not set.")
    }
    else {
      w@aggregated <- filterMultiplicity(w, archivename, 
                                         mode, multiplicityFilter = settings$multiplicityFilter)
      w@aggregated <- processProblematicPeaks(w, mode, 
                                              archivename)
      if (!is.na(archivename)) 
        archiveResults(w, paste(archivename, "_RF.RData", 
                                sep = ""), settings)
    }
  }
  message("msmsWorkflow: Done.")
  return(w)
}

msmsRead1 <- function(w, filetable = NULL, files = NULL, cpdids = NULL, 
                      readMethod, mode, confirmMode = FALSE, useRtLimit = TRUE, 
                      Args = NULL, settings = getOption("RMassBank"),
                      progressbar = "progressBarHook", MSe = FALSE, plots = FALSE){
  #.checkMbSettings()
  ##Read the files and cpdids according to the definition
  ##All cases are silently accepted, as long as they can be handled according to one definition
  if(!any(mode %in% c("pH","pNa","pM","pNH4","mH","mFA","mM",""))) stop(paste("The ionization mode", mode, "is unknown."))
  
  if(is.null(filetable)){
    ##If no filetable is supplied, filenames must be named explicitly
    if(is.null(files))
      stop("Please supply the files")
    
    ##Assign the filenames to the workspace
    w@files <- unlist(files)
    
    ##If no filetable is supplied, cpdids must be delivered explicitly or implicitly within the filenames
    if(is.null(cpdids)){
      splitfn <- strsplit(files,"_")
      splitsfn <- sapply(splitfn, function(x) x[length(x)-1])
      if(suppressWarnings(any(is.na(as.numeric(splitsfn)[1]))))
        stop("Please supply the cpdids corresponding to the files in the filetable or the filenames")
      cpdids <- splitsfn
    }
  } else{
    ##If a filetable is supplied read it
    tab <- read.csv(filetable, stringsAsFactors = FALSE)
    w@files <- tab[,"Files"]
    cpdids <- tab[,"ID"]
  }
  
  ##If there's more cpdids than filenames or the other way around, then abort
  if(length(w@files) != length(cpdids)){
    stop("There are a different number of cpdids than files")
  }
  
  if(!(readMethod %in% c("mzR","peaklist","xcms","minimal"))){
    stop("The supplied method does not exist")
  }
  
  if(!all(file.exists(w@files))){
    stop("The supplied files ", paste(w@files[!file.exists(w@files)]), " don't exist")
  }
  
  # na.ids <- which(is.na(sapply(cpdids, findSmiles)))
  
  # if(length(na.ids)){
  # stop("The supplied compound ids ", paste(cpdids[na.ids], collapse=" "), " don't have a corresponding smiles entry. Maybe they are missing from the compound list")
  # }
  
  ##This should work
  if(readMethod == "minimal"){
    ##Edit options
    opt <- getOption("RMassBank")
    opt$recalibrator$MS1 <- "recalibrate.identity"
    opt$recalibrator$MS2 <- "recalibrate.identity"
    opt$add_annotation==FALSE
    options(RMassBank=opt)
    ##Edit analyzemethod
    analyzeMethod <- "intensity"
  }
  
  if(readMethod == "mzR"){
    ##Progressbar
    nLen <- length(w@files)
    nProg <- 0
    pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
    
    count <- 1
    envir <- environment()
    w@spectra <- as(lapply(w@files, function(fileName) {
      
      # Find compound ID
      cpdID <- cpdids[count]
      retrieval <- findLevel(cpdID,TRUE)
      # Set counter up
      envir$count <- envir$count + 1
      
      # Retrieve spectrum data
      
      spec <- findMsMsHR1(fileName = fileName, 
                          cpdID = cpdID, mode = mode, confirmMode = confirmMode, useRtLimit = useRtLimit,
                          ppmFine = settings$findMsMsRawSettings$ppmFine,
                          mzCoarse = settings$findMsMsRawSettings$mzCoarse,
                          fillPrecursorScan = settings$findMsMsRawSettings$fillPrecursorScan,
                          rtMargin = settings$rtMargin,
                          deprofile = settings$deprofile, retrieval=retrieval)
      
      gc()
      
      # Progress:
      nProg <<- nProg + 1
      pb <- do.call(progressbar, list(object=pb, value= nProg))
      
      return(spec)
    } ), "SimpleList")
    names(w@spectra) <- basename(as.character(w@files))
    return(w)
  }
  
  ##xcms-readmethod 
  if(readMethod == "xcms"){
    
    ##Load libraries
    requireNamespace("xcms",quietly=TRUE)
    requireNamespace("CAMERA",quietly=TRUE)
    
    ##Find unique files and cpdIDs
    ufiles <- unique(w@files)
    uIDs <- unique(cpdids)
    nLen <- length(ufiles)
    
    ##Progressbar
    nProg <- 0
    pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
    i <- 1
    
    ##Routine for the case of multiple cpdIDs per file
    if(length(uIDs) > length(ufiles)){
      w@spectra <- as(unlist(lapply(ufiles, function(currentFile){
        fileIDs <- cpdids[which(w@files == currentFile)]
        spec <- findMsMsHRperxcms(currentFile, fileIDs, mode=mode, findPeaksArgs=Args, plots, MSe = MSe)
        gc()
        
        # Progress:
        nProg <<- nProg + 1
        pb <- do.call(progressbar, list(object=pb, value= nProg))
        
        return(spec)
      }),FALSE),"SimpleList")
      return(w)
    }
    
    ##Routine for the other cases
    w@spectra <- as(lapply(uIDs, function(ID){
      # Find files corresponding to the compoundID
      currentFile <- w@files[which(cpdids == ID)]
      
      # Retrieve spectrum data
      spec <- findMsMsHRperxcms(currentFile, ID, mode=mode, findPeaksArgs=Args, plots, MSe = MSe)
      gc()
      
      # Progress:
      nProg <<- nProg + 1
      pb <- do.call(progressbar, list(object=pb, value= nProg))
      
      return(spec)
    }),"SimpleList")
    ##If there are more files than unique cpdIDs, only remember the first file for every cpdID
    w@files <- w@files[sapply(uIDs, function(ID){
      return(which(cpdids == ID)[1])
    })]
    return(w)
  }
  
  ##Peaklist-readmethod 
  if((readMethod == "peaklist") || (readMethod=="minimal")){
    w <- createSpecsFromPeaklists(w, cpdids, filenames=w@files, mode=mode)
    uIDs <- unique(cpdids)
    files <- list()
    
    for(i in 1:length(uIDs)){
      indices <- sapply(cpdids,function(a){return(uIDs[i] %in% a)})
      files[[i]] <- w@files[indices]
    }
    
    w@files <- sapply(files,function(file){return(file[1])})
    message("Peaks read")
    return(w)
  }
}

findMsMsHR1 <- function(fileName = NULL, msRaw = NULL, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE,
                        ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
                        mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
                        fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
                        rtMargin = getOption("RMassBank")$rtMargin,
                        deprofile = getOption("RMassBank")$deprofile,
                        headerCache = NULL,
                        peaksCache = NULL,
                        retrieval="standard"){
  
  # access data directly for finding the MS/MS data. This is done using
  # mzR.
  if(!is.null(fileName) & !is.null(msRaw))
    stop("Both MS raw data and MS filename given. Only one can be handled at the same time.")
  if(!is.null(fileName))
    msRaw <- openMSfile(fileName)
  
  mzLimits <- findMz(cpdID, mode, retrieval=retrieval)
  mz <- mzLimits$mzCenter
  limit.fine <- ppm(mz, ppmFine, p=TRUE)
  if(!useRtLimit)
    rtLimits <- NA
  else
  {
    dbRt <- findRt(cpdID)
    rtLimits <- c(dbRt$RT - rtMargin, dbRt$RT + rtMargin) * 60
  }
  spectra <- findMsMsHR.mass1(msRaw, mz, mzCoarse, limit.fine, rtLimits, confirmMode + 1,headerCache
                              ,fillPrecursorScan, deprofile, peaksCache, cpdID)
  
  # check whether a) spectrum was found and b) enough spectra were found
  if(length(spectra) < (confirmMode + 1))
    sp <- new("RmbSpectraSet", found=FALSE)
  else
    sp <- spectra[[confirmMode + 1]]
  
  #sp@mz <- mzLimits
  sp@id <- as.character(as.integer(cpdID))
  #sp@name <- findName(cpdID)
  ENV <- environment()
  if(retrieval == "unknown"){
    sp@formula <- ""
  } else{
    sp@formula <- findFormula(cpdID, retrieval=retrieval)
  }
  sp@mode <- mode
  
  # If we had to open the file, we have to close it again
  if(!is.null(fileName))
    close(msRaw)
  
  return(sp)
}

findMsMsHR.mass1 <- function(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
                             headerCache = NULL, fillPrecursorScan = FALSE,
                             deprofile = getOption("RMassBank")$deprofile, peaksCache = NULL, cpdID = NA)
{
  eic <- findEIC(msRaw, mz, limit.fine, rtLimits, headerCache=headerCache, 
                 peaksCache=peaksCache)

  #	if(!is.na(rtLimits))
  #	{  
  #		eic <- subset(eic, rt >= rtLimits[[1]] & rt <= rtLimits[[2]])
  #	}
  if(!is.null(headerCache))
    headerData <- headerCache
  else
    headerData <- as.data.frame(header(msRaw))
  
  ###If no precursor scan number, fill the number
  if(length(unique(headerData$precursorScanNum)) == 1){
    fillPrecursorScan <- TRUE
  }
  
  if(fillPrecursorScan == TRUE)
  {
    # reset the precursor scan number. first set to NA, then
    # carry forward the precursor scan number from the last parent scan
    headerData$precursorScanNum <- NA
    headerData[which(headerData$msLevel == 1),"precursorScanNum"] <-
      headerData[which(headerData$msLevel == 1),"acquisitionNum"]
    headerData[,"precursorScanNum"] <- .locf(headerData[,"precursorScanNum"])
    # Clear the actual MS1 precursor scan number again
    headerData[which(headerData$msLevel == 1),"precursorScanNum"] <- 0
  }
  # bugfix 201803: PRM scans that were performed before the first full scan (found in some files)
  headerData <- headerData[
    !((headerData$msLevel == 2) & (headerData$precursorScanNum == 0)),,drop=FALSE
    ]
  # Find MS2 spectra with precursors which are in the allowed 
  # scan filter (coarse limit) range

  findValidPrecursors <- headerData[
    (headerData$precursorMZ > mz - limit.coarse) &
      (headerData$precursorMZ < mz + limit.coarse),,drop=FALSE]
  # Find the precursors for the found spectra
  validPrecursors <- unique(findValidPrecursors$precursorScanNum)
  # check whether the precursors are real: must be within fine limits!
  # previously even "bad" precursors were taken. e.g. 1-benzylpiperazine
  validPrecursors = validPrecursors[!is.na(validPrecursors)]

  which_OK <- lapply(validPrecursors, function(pscan)
  {
    pplist <- as.data.frame(
      peaks(msRaw, which(headerData$acquisitionNum == pscan)))
    colnames(pplist) <- c("mz","int")
    pplist <- pplist[(pplist$mz >= mz -limit.fine)
                     & (pplist$mz <= mz + limit.fine),,drop=FALSE]
    if(nrow(pplist) > 0)
      return(TRUE)
    return(FALSE)
  })
  validPrecursors <- validPrecursors[which(which_OK==TRUE)]
  if(length(validPrecursors) == 0){
    if(!is.na(cpdID))
      warning(paste0("No precursor was detected for compound, ", cpdID, " with m/z ", mz, ". Please check the mass and retention time window."))
    else
      warning(paste0("No precursor was detected for m/z ", mz, ". Please check the mass and retention time window."))
  }
  # Crop the "EIC" to the valid precursor scans
  eic <- eic[eic$scan %in% validPrecursors,]
  # Order by intensity, descending
  eic <- eic[order(eic$intensity, decreasing=TRUE),]
  
  if(nrow(eic) == 0)
    return(list(
      new("RmbSpectraSet",
          found=FALSE)))
  if(!is.na(maxCount))
  {
    spectraCount <- min(maxCount, nrow(eic))
    eic <- eic[1:spectraCount,,drop=FALSE]
  }

  # Construct all spectra groups in decreasing intensity order
  spectra <- lapply(eic$scan, function(masterScan)
  {
    masterHeader <- headerData[headerData$acquisitionNum == masterScan,,drop=FALSE]
    
    childHeaders <- headerData[which((headerData$precursorScanNum == masterScan) 
                               & (headerData$precursorMZ > mz - limit.coarse) 
                               & (headerData$precursorMZ < mz + limit.coarse)),,drop=FALSE]
    
    # Fix 9.10.17: headers now include non-numeric columns, leading to errors in data conversion.
    # Remove non-numeric columns
    headerCols <- colnames(masterHeader)
    headerCols <- headerCols[unlist(lapply(headerCols, function(col) is.numeric(masterHeader[,col])))]
    masterHeader <- masterHeader[,headerCols,drop=FALSE]
    childHeaders <- childHeaders[,headerCols,drop=FALSE]
    
    childScans <- childHeaders$seqNum
    
    msPeaks <- peaks(msRaw, masterHeader$seqNum)

    # if deprofile option is set: run deprofiling
    deprofile.setting <- deprofile
    if(!is.na(deprofile.setting))
      msPeaks <- deprofile.scan(
        msPeaks, method = deprofile.setting, noise = NA, colnames = FALSE
      )
    colnames(msPeaks) <- c("mz","int")

    msmsSpecs <- apply(childHeaders, 1, function(line)
    {
      pks <- peaks(msRaw, line["seqNum"])
  
      if(!is.na(deprofile.setting))
      {								
        pks <- deprofile.scan(
          pks, method = deprofile.setting, noise = NA, colnames = FALSE
        )
      }
      
      new("RmbSpectrum2",
          mz = pks[,1],
          intensity = pks[,2],
          precScanNum = as.integer(line["precursorScanNum"]),
          precursorMz = line["precursorMZ"],
          precursorIntensity = line["precursorIntensity"],
          precursorCharge = as.integer(line["precursorCharge"]),
          collisionEnergy = line["collisionEnergy"],
          tic = line["totIonCurrent"],
          peaksCount = line["peaksCount"],
          rt = line["retentionTime"],
          acquisitionNum = as.integer(line["seqNum"]),
          centroided = TRUE
      )
    })
    msmsSpecs <- as(do.call(c, msmsSpecs), "SimpleList")
    
    # build the new objects
    masterSpec <- new("Spectrum1",
                      mz = msPeaks[,"mz"],
                      intensity = msPeaks[,"int"],
                      polarity = as.integer(masterHeader$polarity),
                      peaksCount = as.integer(masterHeader$peaksCount),
                      rt = masterHeader$retentionTime,
                      acquisitionNum = as.integer(masterHeader$seqNum),
                      tic = masterHeader$totIonCurrent,
                      centroided = TRUE
    )
    
    spectraSet <- new("RmbSpectraSet",
                      parent = masterSpec,
                      children = msmsSpecs,
                      found = TRUE,
                      #complete = NA,
                      #empty = NA,
                      #formula = character(),
                      mz = mz
                      #name = character(),
                      #annotations = list()
    )
    return(spectraSet)
  })
  names(spectra) <- eic$acquisitionNum
  return(spectra)
}


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
    
    filter = which(sp1[,2]>=min_relative & sp1[,1]<mz0-1)
    if (normalized){sp = sp1}  
    sp = sp[filter,,drop=FALSE]
    
    # Check validity:
    
    if (nrow(sp)>0 & checked){
      sp = sp[order(sp[,1]),,drop=FALSE]
      denoised_spectrum = sp
    }
  }
  return(denoised_spectrum)
}

denoise_ms1_spectrum<-function(sp, mz0, max_peak, min_relative, normalized = T){
  
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
    
    filter = which(sp1[,2]>=min_relative & sp1[,1]>=mz0-0.5 & sp1[,1]<=mz0+10)
    if (normalized){sp = sp1}  
    sp = sp[filter,,drop=FALSE]
    
    # Check validity:
    
    if (nrow(sp)>=2 & checked){
      sp = sp[order(sp[,1]),]
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


