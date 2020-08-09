#' Read and combine targeted MS1 and MS2 scans from one LC-MS/MS file using compMS2Miner algorithm
#'
#' Function used by library_generator
#' 
#' @importFrom BiocGenerics basename
#' @importFrom utils flush.console txtProgressBar setTxtProgressBar
#' @importFrom tcltk tk_choose.dir tclvalue tkgetOpenFile
#' @importFrom mzR openMSfile header peaks
#' 
#' @export
#' 
process_compMS2Miner<-function(mzdatafiles = NULL, ref = NULL, polarity = c("Positive", "Negative"), include.MS1 = FALSE,
                               rt_search = 10, rt_gap = 30, ppm_search = 10, mz_search = 0.01,
                               baseline= 1000, relative = 5, max_peaks = 200, normalized=T){
  
  options(stringsAsFactors = F)
  options(warn=-1)
  library(compMS2Miner)
  
  # importFrom compMS2Miner deconvNoise combineMS2 compSpectra metaData filePaths
  
  output_metadata = c()
  output_sp = list()
  
  ###########################
  ######Prepare metadata#####
  ###########################
  
  if ("FILENAME" %in% colnames(ref)){
    valid = c(which(basename(ref$FILENAME) == basename(mzdatafiles)),  which(ref$FILENAME =="N/A"))
    ref = ref[valid,,drop=FALSE]
  }
  
  if (nrow(ref)>0){
    
    ID0 = ref$ID #  Backup ID
  
    peakTable = cbind.data.frame(EICno = 1:nrow(ref), mzmed = ref$PEPMASS, rtmed = ref$RT*60, abundance = 100)

    if (polarity=="Positive"){polarity = "pos"}
    if (polarity=="Negative"){polarity = "neg"}
  
    #####################
    ######Processing#####
    #####################
    
    save(peakTable, mzdatafiles, polarity, ppm_search, rt_search, baseline, file = "temp.RData")
    
    compMS2Demo <- compMS2Construct1(MS1features = peakTable, MS2files = mzdatafiles, mode = polarity, 
                                   precursorPpm = ppm_search, ret = rt_search, TICfilter = baseline)
  
    if (length(metaData(compMS2Demo))>0){compMS2Demo <- deconvNoise(compMS2Demo, "DNF")}
    if (length(metaData(compMS2Demo))>0){compMS2Demo <- combineMS2(compMS2Demo, "Ions")}
    if (length(metaData(compMS2Demo))>0){compMS2Demo <- combineMS2(compMS2Demo, "Spectra", specSimFilter=0.5)}
    if (length(metaData(compMS2Demo))>1){
      compMS2Demo <- combineMS2(compMS2Demo, 'removeContam', maxRtGap=rt_gap, minSimScore= 0.5, ms1Abs=mz_search)}

    NFeatures = length(metaData(compMS2Demo))
    IndFeatures = as.numeric(gsub("CC_", "", names(compMS2Demo@compSpectra)))
    IDFeatures =  ID0[IndFeatures]

    #########################
    ######Transformation ####
    #########################

    if (NFeatures>0){
  
      for (i in 1:NFeatures){
    
      temp_sp = list()
      
      old_metadata = ref[IndFeatures[i],,drop=FALSE]
      to_remove = match(c("PEPMASS", "RT"), colnames(ref))
      old_metadata = old_metadata[, -to_remove] # Remove PEPASS and RT
      new_metadata = data.frame(compMS2Demo@metaData[[i]])
      sp0 = data.matrix(compMS2Demo@compSpectra[[i]])
    
      ind_m0 = grep("_precursorMz", colnames(new_metadata))
      ind_rt = grep("_retentionTime", colnames(new_metadata))
      ind_tic = grep("_TIC", colnames(new_metadata))[1]
      ind_dev = grep("_ppmDiff", colnames(new_metadata))
      ind_ms2_scan_nb = grep("_acquisitionNum", colnames(new_metadata))

      # MS2 Metadata:
      
      new_metadata_ms2 = new_metadata[, c(ind_m0, ind_rt, ind_tic, ind_dev, ind_ms2_scan_nb)]
      colnames(new_metadata_ms2) = c("PEPMASS", "RT", "TIC", "PEPMASS_DEV", "SCAN_NUMBER")
      new_metadata_ms2 = new_metadata_ms2[which.max(new_metadata_ms2$TIC),,drop=FALSE] # Save scan with highest TIC
       
      new_metadata_ms2$PEPMASS = round(new_metadata_ms2$PEPMASS, 5)
      new_metadata_ms2$RT = round(mean(new_metadata_ms2$RT)/60, 2)
      new_metadata_ms2$PEPMASS_DEV = round(abs(new_metadata_ms2$PEPMASS_DEV), 2)
      
      temp_metadata = cbind.data.frame(new_metadata_ms2[,c(1:2)], old_metadata, 
                        FILENAME = basename(mzdatafiles), MSLEVEL = 2, new_metadata_ms2[,3:5])
      
      # MS2 Spectrum:
      
      colnames(sp0) = NULL
      temp_sp[[1]] = denoise_ms2_spectrum(sp0, ref$PEPMASS[IndFeatures[i]], max_peaks, relative, normalized)

      if (include.MS1){

        ind_ms1_scan_nb = grep("_precursorScanNum", colnames(new_metadata))
        ind_ms1_tic = grep("_precursorIntensity", colnames(new_metadata))
        ind_ms1 = grep("_isoMass", colnames(new_metadata))
        ind_ms1_int = grep("_isoInt", colnames(new_metadata))
        
        # MS1 Metadata:
        
        selected_metadata = new_metadata[which.max(new_metadata[,ind_ms1_tic]),, drop=FALSE]
             
        new_metadata_ms1 = temp_metadata
        new_metadata_ms1$MSLEVEL = 1
        new_metadata_ms1$TIC = selected_metadata[1,ind_ms1_tic]
        new_metadata_ms1$SCAN_NUMBER = selected_metadata[1,ind_ms1_scan_nb]

        temp_metadata = rbind.data.frame(temp_metadata, new_metadata_ms1)
        
        # MS1 Spectrum:
        
        temp_mass_peaks = selected_metadata[1, ind_ms1]
        temp_mass_peaks = as.numeric(strsplit(temp_mass_peaks, ";")[[1]])
        temp_mass_intensity = selected_metadata[1, ind_ms1_int]
        temp_mass_intensity = as.numeric(strsplit(temp_mass_intensity, ";")[[1]])
        temp_sp_ms1 = data.matrix(cbind(temp_mass_peaks, temp_mass_intensity))
        colnames(temp_sp_ms1) =NULL
        temp_sp_ms1 = denoise_ms1_spectrum(temp_sp_ms1, ref$PEPMASS[IndFeatures[i]], max_peaks, relative, normalized)
        temp_sp[[2]] = temp_sp_ms1
      }
    
      output_metadata = rbind.data.frame(output_metadata, temp_metadata)
      output_sp = c(output_sp, temp_sp)
  }
    
  ###########################
  ####Filter empty spectra###
  ###########################
      
  filter = which(sapply(output_sp, nrow)>0)
  output_metadata = output_metadata[filter,,drop=FALSE]
  output_sp = output_sp[filter]
 }}  
  
  return(list(sp=output_sp,metadata=output_metadata))
}

############################
### Internal functions:####
###########################

compMS2Construct1 <-  function(MS1features = NULL, msDataDir = NULL, MS2files=NULL, 
                               mode = "pos", precursorPpm = 10, ret = 10, 
                               TICfilter = 10000, minPeaks=1, isoWid=4, 
                               verbose=TRUE){
  
  options(stringsAsFactors = F)
  options(warn=-1)
  
  message("creating compMS2 object in ", ifelse(mode == "pos","positive", "negative"), 
          " ionisation mode")
  flush.console()
  
  # new compMS2 object 
  object <- new("compMS2")
  # set global options
  # if is.null msDataDir select .mzXML/.mzML/.mgf file containing directory
  if(is.null(MS2files)){
    if(is.null(msDataDir)){
      message("Select your .mzXML/.mzML/.mgf data directory")
      flush.console()
      msDataDir <- tk_choose.dir(default = "",  
                                        caption = "1. Select your .mzXML/.mzML/.mgf data directory")
    }
    
    # identify all mzXML/.mzML files in raw-data directory
    MS2files <- list.files(path=msDataDir, pattern = "\\.mzXML$|\\.mzML$|\\.mgf$", 
                           full.names=TRUE)
  }
  fileTypeTmp <- unique(gsub('.+\\.', '', basename(MS2files)))
  if(length(fileTypeTmp) > 1){
    stop('only 1 file type is permitted in the MS data directory (msDataDir):\n',
         'The following file types were found:\n', paste0('.', fileTypeTmp, '\n'))
  }
  message(length(MS2files), 
          paste0(" MS2 (.", fileTypeTmp, ") file(s) were detected within the directory..."))
  flush.console()
  adducts <- FALSE
  # if file type is mgf
  if(fileTypeTmp == 'mgf'){
    # file names
    fileNames <- gsub('\\.mgf$', '', basename(MS2files))
    # remove hypens
    fileNames <- gsub('-', '_', fileNames)
    # spectra
    compSpectra(object) <- vector('list', length(MS2files))
    names(compSpectra(object)) <- fileNames
    # metaData
    metaData(object) <- vector('list', length(MS2files))
    names(metaData(object)) <- fileNames
    message('reading mgf files...\n')
    flush.console()
    if(verbose == TRUE){ pb <- txtProgressBar(max=length(MS2files), style = 3)}
    
    for(i in 1:length(MS2files)){
      if(verbose==TRUE){setTxtProgressBar(pb, i)}
      mgfFile <- readLines(MS2files[i])
      # spectrum
      spectrumTmp <- mgfFile[grep('^[0-9]', mgfFile)]
      spectrumTmp <- do.call(rbind, strsplit(spectrumTmp, '\\t'))
      if(nrow(spectrumTmp) > 1){
        spectrumTmp <- data.frame(apply(spectrumTmp, 2, as.numeric), 
                                  stringsAsFactors=FALSE)
      } else {
        spectrumTmp <- as.numeric(spectrumTmp[1, ])
        spectrumTmp <- data.frame(spectrumTmp[1], spectrumTmp[2], stringsAsFactors=FALSE)
      }
      colnames(spectrumTmp) <- c('mass', 'intensity')
      compSpectra(object)[[i]] <- spectrumTmp
      # metaData
      massIndx <- grep('PEPMASS=', mgfFile) 
      rtIndx <- grep('RTINSECONDS=', mgfFile) 
      metaDataTmp <- list(as.numeric(gsub('PEPMASS=', '', mgfFile[massIndx])), 
                          as.numeric(gsub('RTINSECONDS=', '', mgfFile[rtIndx])),
                          0)
      names(metaDataTmp) <- paste0(fileNames[i], c('_MS1_mz', '_MS1_RT', '_precursorIntensity'))
      metaData(object)[[i]] <- metaDataTmp
    }
  } else { # if .mzXML or .mzML
    # if is.null MS1features select MS1 feature table file
    if(is.null(MS1features)){
      message("Select your MS1 feature (.csv) file")
      flush.console()
      
      MS1features <- tclvalue(tkgetOpenFile(
        initialdir = msDataDir,
        filetypes = "{{Comma delimited} {.csv}} {{All files} *}", 
        title = "2. select the MS1 features (.csv) file you wish to match"))
    }
    
    if(is.character(MS1features)){
      # read in MS1features
      message("Reading MS1 feature table...")
      flush.console()
      MS1features <- as.data.frame(data.table::fread(MS1features, sep=",", 
                                                     header=TRUE, stringsAsFactors=FALSE))
      message("...Done")
      flush.console()
    }
    
    if(!is.data.frame(MS1features)){
      stop('The MS1features object is not a data.frame.')
    }
    # error handling
    if(!is.integer(MS1features[, 1])){
      stop('The first column of the MS1 feature table must be an integer (EIC/unique number).')
    }
    if(!is.numeric(MS1features[, 2])){
      stop('The second column of the MS1 feature table must be a numeric (mass-to-charge ratio).')
    }
    if(!is.numeric(MS1features[, 3])){
      stop('The third column of the MS1 feature table must be a numeric (retention time in seconds).')
    }
    # sort dataframe by unique ID/ EIC number
    MS1features <- MS1features[order(MS1features[, 1]), ]
    row.names(MS1features) <- seq(1, nrow(MS1features), 1)
    # see if the 4 th column contains adduct information
    adducts <- any(grepl('M\\+|M\\-', MS1features[, 4]))
    if(adducts == TRUE){
      message('Adducts/fragments were detected in the 4th column of the MS1features table. These may be used for subsequent stages of the workflow.\n')
      flush.console()
    }
    
    # create list to store results
    Results <- vector("list", length(MS2files))
    
    for(i in 1:length(MS2files)){
      Res.tmp <- compMS2Create1(MS2file=MS2files[i], TICfilter=TICfilter, 
                                MS1features=MS1features, precursorPpm=precursorPpm,
                                ret=ret, adducts=adducts, isoWid=isoWid)

      Results[[i]] <- Res.tmp
    }

    Results <- unlist(Results, recursive = FALSE)

    if (length(Results)>0){
      compSpectra(object) <- lapply(Results, function(x) x$spectra)
      metaData(object) <- lapply(Results, function(x) x$metaData)
    } #else {
      #compSpectra(object) <- lapply(Results, function(x) x$spectra)
      #metaData(object) <- lapply(Results, function(x) x$metaData) 
    #}

    # check all compSpectra are not empty
    indxTmp <- sapply(compSpectra(object), nrow) > minPeaks

    if(any(indxTmp == FALSE)){
      compSpectra(object) <- compSpectra(object)[indxTmp]
      metaData(object) <- metaData(object)[indxTmp]
    }
  }
  
  # add file paths and parameters
  filePaths(object) <- MS2files
  Parameters(object) <- data.frame(nCores=0,
                                   mode=mode, precursorPpm=precursorPpm,
                                   ret=ret, TICfilter=TICfilter, 
                                   fileType=fileTypeTmp, adducts=adducts, 
                                   isoWid=isoWid, stringsAsFactors=FALSE)
  return(object)
} # end CompMS2obj function

compMS2Create1 <- function(MS2file = NULL, MS1features = NULL, 
                          TICfilter = 10000, precursorPpm = 10, ret = 10, 
                          adducts=FALSE, isoWid=4){
  
  # MS2 file name
  MS2fileName <- basename(MS2file)
  
  message(paste0("Reading ", MS2fileName, "..."))
  flush.console()
  # read MS2 file
  MS2file <- openMSfile(MS2file)
  
  message("...DONE")
  flush.console()
  
  message("extracting metaData from MS2 file")
  flush.console()
  
  metaData <- header(MS2file)
  metaData <- metaData[, c('msLevel', 'precursorMZ', 'retentionTime',
                           'totIonCurrent', 'precursorIntensity',
                           'collisionEnergy', 'basePeakMZ', 'basePeakIntensity',
                           'acquisitionNum', 'precursorScanNum')]
  metaData$TICaboveFilter <- {metaData$totIonCurrent >= TICfilter} * 1
  colnames(metaData) <- c("MS.scanType", "precursorMz", "retentionTime", "TIC", 
                          "precursorIntensity", 
                          "collisionEnergy", "basePeakMz", "basePeakIntensity",
                          "acquisitionNum","precursorScanNum", "TICaboveFilter")
  
  metaData <- metaData[, c("MS.scanType", "precursorMz", "retentionTime", "TIC", 
                           "TICaboveFilter", "precursorIntensity", 
                           "collisionEnergy", "basePeakMz", "basePeakIntensity",
                           "acquisitionNum","precursorScanNum")]
  
  message("...DONE")
  flush.console()
  
  # cond if no MS2 level scans detected
  if(all(metaData$MS.scanType == 1)){
    # no MS2 level scans detected cond
    warning(paste0("No MS2 levels scans within ", MS2fileName, ",  check that the 
                   file has been converted to the mzXML format correctly."), 
            immediate.=TRUE)
    flush.console()
    message("...moving to next MS2 file")
    flush.console()
  } else {
    # remove all MS/MS scans where the TIC is less than the minimum TIC threshold 
    # set by the user
    message(paste0("Of a total of ", length(which(metaData$MS.scanType == 2)),
                   " MS2 spectra..."))
    flush.console()
    # index ms2 scan and above TIC filter
    metaData$MS2TICfilt.indx <- (metaData$MS.scanType == 2 & 
                                   metaData$TICaboveFilter == 1) * 1
    nAboveTIC <- length(which(metaData$MS2TICfilt.indx == 1))
    message(paste0(nAboveTIC, " MS2 spectra were above the TIC filter of ", 
                   TICfilter))
    flush.console()
    # cond if no scan above the TIC filter
    if(length(nAboveTIC) == 0){ 
      warning(paste0("No MS2 levels scans above TIC filter of ", TICfilter, " in ", 
                     MS2fileName, ",  reduce the TIC filter parameter or check that 
                     the file has been converted to the mzXML format correctly."), 
              immediate.=TRUE)
      flush.console()
      message("...moving to next MS2 file")
      flush.console()
    } else { 
      
      message("matching MS1 peak table to precursors of MS2 spectra...")
      flush.console()
      
      #save(metaData, MS2file, adducts, isoWid, MS1features, file = "tmp.RData")
      # mapply MS1 feature match

      MS1MS2match <- mapply(MS1MatchSpectra1, EIC=MS1features[, 1], 
                            mz=MS1features[, 2], RT=MS1features[, 3], 
                            adduct=MS1features[, 4],
                            precursorPpm=precursorPpm, ret=ret, 
                            MoreArgs=list(metaData=metaData, 
                                          MS2file=MS2file, adducts=adducts, 
                                          isoWid=isoWid))
      # for(i in 1:nrow(MS1features)){
      #   tmp <- MS1MatchSpectra(EIC=MS1features[i, 1], 
      #                          mz=MS1features[i, 2], RT=MS1features[i, 3], 
      #                          adduct=MS1features[i, 4],
      #                          precursorPpm=precursorPpm, ret=ret,
      #                          metaData=metaData, 
      #                          MS2file=MS2file, adducts=adducts, 
      #                          isoWid=isoWid)
      #   # EIC=MS1features[i, 1]; 
      #   # mz=MS1features[i, 2]; RT=MS1features[i, 3]; 
      #   # adduct=MS1features[i, 4];
      # }

      message("...done")
      flush.console()
      
      if (is.data.frame(MS1MS2match[[1]])){
       new_MS1MS2match = list()
       new_MS1MS2match[[1]] = list()
       new_MS1MS2match[[1]]$spectra = MS1MS2match[[1]]
       new_MS1MS2match[[1]]$metaData = MS1MS2match[[2]]
       MS1MS2match = new_MS1MS2match
      }
        
      match.indx <- which(sapply(MS1MS2match, length) == 2)
    
      # calculate composite spectra
      message(paste0(length(match.indx), " MS1 features were matched to MS2 precursors"))
      flush.console()
      names(MS1MS2match) <- paste0(MS2fileName, "_", MS1features[match.indx, 1])

      # check for chimeric spectra and isotopes
      
      MS1MS2match <- MS1MS2match[match.indx]
      
      return(MS1MS2match)
      #Results[names(MS1MS2match)] <- MS1MS2match
      # close(MS2file)
      #time[i] <- (proc.time() - pmt)[["elapsed"]] 
    } # cond if no scan above the TIC filter
  } # cond if no MS2 level scans detected
} # end func


MS1MatchSpectra1 <- function(metaData=NULL, MS2file=NULL, mz=NULL, RT=NULL, 
                            EIC=NULL, adduct=NULL, precursorPpm = 10, ret = 10, 
                            adducts=FALSE, isoWid=4){
  #error handling
  if(is.null(metaData)){
    stop("metaData value is missing with no default")
  } else if(is.null(MS2file)){
    stop("MS2file is missing with no default")
  } else if(is.null(mz)){
    stop("mz value is missing with no default")
  } else if(is.null(RT)){
    stop("RT value is missing with no default")
  } else if(is.null(EIC)){
    stop("EIC value is missing with no default")
  } else {
    # match MS1 mass by ppm tolerance to precursor m/z
    # match MS1 RT in seconds to precursor RT
    MS1.match <- which({metaData$precursorMz < 
        mz+(mz/1E06)*precursorPpm &
        metaData$precursorMz > 
        mz-(mz/1E06)*precursorPpm &
        metaData$retentionTime < RT+ret &
        metaData$retentionTime > RT-ret &
        metaData$MS2TICfilt.indx == 1} == TRUE)
    
    if(length(MS1.match) == 0){
      MS1.match <- 0
      return(MS1.match)
      
    } else  if(length(MS1.match) == 1){
      
      spectra <- data.frame(peaks(MS2file, MS1.match))
      colnames(spectra) <- c("mass", "intensity")
      #       spectra$scanSeqNum <- mzR::header(MS2file, MS1.match)$seqNum
      metaData.df <- metaData[MS1.match, , drop = FALSE]
      
      # ms1 scans
      
      precursorScans <- unique(metaData.df$precursorScanNum)
      
      # voltage ramps
      if(all(precursorScans == 0)){
        precursorScans <- metaData$acquisitionNum[MS1.match - 1]
        metaData.df$precursorScanNum <- precursorScans
      }
      ms1ScanIdx <- match(metaData$acquisitionNum, precursorScans)
      # subset if precursor scan missing
      metaData.df <- metaData.df[ms1ScanIdx[!is.na(ms1ScanIdx)], , drop=FALSE]
      precursorScans <- precursorScans[ms1ScanIdx[!is.na(ms1ScanIdx)]]
      ms1ScanIdx <- which(!is.na(ms1ScanIdx))
      ms1Prec <- unique(metaData.df$precursorMz)[1]
      isoChim <- data.frame(matrix(0, ncol=4, nrow=length(ms1ScanIdx)), stringsAsFactors = FALSE)
      colnames(isoChim) <-  c('isoMass', 'isoInt', 'possChim', 'maxInterIons')
      for(sc in 1:length(ms1ScanIdx)){
        x <- peaks(MS2file, ms1ScanIdx[sc])
        idxTmp <-  x[, 1] < {ms1Prec + 2.5} & x[, 1] > {ms1Prec - 2.5}
        if(sum(idxTmp) > 2){
          ms1SpecTmp <- x[idxTmp, , drop=FALSE]
          prIdx <- which.min(abs(ms1SpecTmp[, 1] - ms1Prec))
          iso1 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[prIdx, 1] + 1}))
          iso2 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[iso1, 1] + 1}))
          isoMass <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 1], 6), collapse = ';')
          isoInt <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 2], 2), collapse = ';')
          idxTmp <- x[, 1] < {ms1Prec + {isoWid/2}} & x[, 1] > {ms1Prec - {isoWid/2}} 
          possChim <- x[idxTmp, , drop=FALSE]
          prIdx <- which.min(abs(possChim[, 1] - ms1Prec))
          iso1 <- which.min(abs(possChim[, 1] - {possChim[prIdx, 1] + 1}))
          iso2 <- which.min(abs(possChim[, 1] - {possChim[iso1, 1] + 1}))
          possChim[, 2] <- {possChim[, 2]/possChim[prIdx, 2]} * 100
          # remove prec and isotopes within 0.5 Da
          idxClose <- possChim[, 1] < {possChim[prIdx, 1] + 0.5} & possChim[, 1] > {possChim[prIdx, 1] - 0.5} | possChim[, 1] < {possChim[iso1, 1] + 0.5} & possChim[, 1] > {possChim[iso1, 1] - 0.5} | possChim[, 1] < {possChim[iso2, 1] + 0.5} & possChim[, 1] > {possChim[iso2, 1] - 0.5}
          possChim <- possChim[idxClose == FALSE, , drop=FALSE]
          if(nrow(possChim) > 0){
            maxInterIons <- max(possChim[, 2])
            possChim <- any(possChim[, 2] >= 50)
          } else {
            maxInterIons <- 0
            possChim <- FALSE 
          }
          isoChim[sc, ] <- c(isoMass, isoInt, possChim, maxInterIons)
        }
      }
      row.names(isoChim) <- precursorScans
      
      metaData.df <- cbind(metaData.df, 
                           isoChim[match(as.character(metaData.df$precursorScanNum),
                                         row.names(isoChim)), ])

      metaData.df$ppmDiff <- ((mz - metaData.df$precursorMz) / mz) * 1E06 
      metaData.df$rtDiff <- metaData.df$retentionTime - RT
      metaData.df$MS1_EICno <- EIC
      metaData.df$MS1_mz <- mz
      metaData.df$MS1_RT <- RT
      
      if(adducts == TRUE){
        metaData.df$MS1_adduct <- adduct
      } else {
        metaData.df$MS1_adduct <- ''  
      }
      MS1.match <- list(spectra = spectra, metaData = metaData.df)    
    } else if(length(MS1.match) > 1){
      metaData.df <- metaData[MS1.match, , drop = FALSE]
      # ms1 scans
      precursorScans <- unique(metaData.df$precursorScanNum)

      if(all(precursorScans == 0)){
        precursorScans <- metaData$acquisitionNum[MS1.match]
        
        #precursorScans <- metaData$acquisitionNum[MS1.match - 1]
        metaData.df$precursorScanNum <- precursorScans
      } 
      ms1ScanIdx <- match(metaData$acquisitionNum, precursorScans)
      if(any(!is.na(ms1ScanIdx))){
        # subset if precursor scan missing
        metaData.df <- metaData.df[ms1ScanIdx[!is.na(ms1ScanIdx)], , drop=FALSE]
        precursorScans <- precursorScans[ms1ScanIdx[!is.na(ms1ScanIdx)]]
        MS1.match <- MS1.match[ms1ScanIdx[!is.na(ms1ScanIdx)]]
        if(length(MS1.match) == 1){
          spectra <- data.frame(peaks(MS2file, MS1.match))
        } else {
          spectra <- data.frame(do.call(rbind, peaks(MS2file, MS1.match)), stringsAsFactors = FALSE)
        }
        colnames(spectra) <- c("mass", "intensity")
        ms1ScanIdx <- which(!is.na(ms1ScanIdx))
        
        ms1Prec <- unique(metaData.df$precursorMz)[1]
        isoChim <- data.frame(matrix(NA, ncol=4, nrow=length(ms1ScanIdx)), stringsAsFactors = FALSE)
        colnames(isoChim) <-  c('isoMass', 'isoInt', 'possChim', 'maxInterIons')
        
        for(sc in 1:length(ms1ScanIdx)){
          x <- peaks(MS2file, ms1ScanIdx[sc])
          idxTmp <-  x[, 1] < {ms1Prec + 2.5} & x[, 1] > {ms1Prec - 2.5}
          if(sum(idxTmp) > 2){
            ms1SpecTmp <- x[idxTmp, , drop=FALSE]
            prIdx <- which.min(abs(ms1SpecTmp[, 1] - ms1Prec))
            iso1 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[prIdx, 1] + 1}))
            iso2 <- which.min(abs(ms1SpecTmp[, 1] - {ms1SpecTmp[iso1, 1] + 1}))
            isoMass <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 1], 6), collapse = ';')
            isoInt <- paste0(round(ms1SpecTmp[c(prIdx, iso1, iso2), 2], 2), collapse = ';')
            idxTmp <- x[, 1] < {ms1Prec + {isoWid/2}} & x[, 1] > {ms1Prec - {isoWid/2}} 
            possChim <- x[idxTmp, , drop=FALSE]
            prIdx <- which.min(abs(possChim[, 1] - ms1Prec))
            iso1 <- which.min(abs(possChim[, 1] - {possChim[prIdx, 1] + 1}))
            iso2 <- which.min(abs(possChim[, 1] - {possChim[iso1, 1] + 1}))
            possChim[, 2] <- {possChim[, 2]/possChim[prIdx, 2]} * 100
            # remove prec and isotopes
            idxClose <- possChim[, 1] < {possChim[prIdx, 1] + 0.5} & possChim[, 1] > {possChim[prIdx, 1] - 0.5} | possChim[, 1] < {possChim[iso1, 1] + 0.5} & possChim[, 1] > {possChim[iso1, 1] - 0.5} | possChim[, 1] < {possChim[iso2, 1] + 0.5} & possChim[, 1] > {possChim[iso2, 1] - 0.5}
            possChim <- possChim[idxClose == FALSE, , drop=FALSE]
            if(nrow(possChim) > 0){
              maxInterIons <- max(possChim[, 2])
              possChim <- any(possChim[, 2] >= 50)
            } else {
              maxInterIons <- 0
              possChim <- FALSE 
            }
            isoChim[sc, ] <- c(isoMass, isoInt, possChim, maxInterIons)
          }
        }
        row.names(isoChim) <- precursorScans
        
        metaData.df <- cbind.data.frame(metaData.df, 
                             isoChim[match(as.numeric(metaData.df$precursorScanNum),
                                           row.names(isoChim)), ])

        metaData.df$ppmDiff <- ((mz - metaData.df$precursorMz) / mz) * 1E06 
        metaData.df$rtDiff <- metaData.df$retentionTime - RT
        metaData.df$MS1_EICno <- EIC
        metaData.df$MS1_mz <- mz
        metaData.df$MS1_RT <- RT
        
        if(adducts == TRUE){
          metaData.df$MS1_adduct <- adduct
        } else {
          metaData.df$MS1_adduct <- ''  
        }
        MS1.match <- list(spectra = spectra, metaData = metaData.df)    
      } else {
        MS1.match <- 0
      }
    }
    return(MS1.match)
  }
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