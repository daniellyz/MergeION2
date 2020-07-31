#' Generating spectral library from raw LC-MS/MS chromatograms
#'
#' The function proposes three data processing algorithms to pick up MS1/MS2 scans from DDA or targeted mode LC-MS/MS data and merge them into a spectral library (new or existing).
#' 
#' @param input_library Character or library object. If character, name of the library into which new scans are added, the file extension must be mgf, msp or RData; please set to empty string "" or NULL if the new library has no dependency with previous ones.
#' @param raw_data_files A character vector of LC-MS/MS file names from which scans are extracted. All files must have be in centroid-mode with mzML, mzXML or cdf extension!
#' @param metadata_file A single character. It should be the metadata file name. The file should be tab, comma or semi-colon separated txt, dat or csv format. For all algorithms, the metadata must contain the column "ID" - a unique structure identifier. The column PEPMASS (targeted precursor mass) must be provided for Default and compMS2Miner. The column RT (targeted retention time in min) must be provided for compMS2Miner and optional for MergeION and RMassBank. Please include the column SMILES (structure identifier) for RMassBank algorithm. If RMassBank is used, the column FILENAME (chromatogram file with mzML, mzXML or cdf extension) must be provided for each compound telling the algorithm from which file compound can be found. Column FILENAME is optional for Default and compMS2Miner. Column ADDUCT is optional for all algorithms, if not provided, all input will be considered as M+H or M-H depending on polarity. Please specify the adduct type if metadata contains both positive and negative ions.
#' @param polarity A single character. Either "Positive" or "Negative". Ion mode of LC-MS/MS files. 
#' @param mslevel A numeric vector. Must contain 2 (if only MS2 scans are extracted) and can be c(1,2) if isotopic pattern in MS1 scans are also extracted. Note: High-quality isotopic patterns in MS1 scans are useful for determining precursor formula!
#' @param add.adduct Logical. If TRUE, additional adduct types will be calculated based on precursor masses of "M+H" and "M-H" adducts in the input metadata: "M+2H", "M+Na","M+K","M+NH4","M+" will be searched for positive ion mode, "M+COO-", "M+Cl" and "M+CH3COO-" for negative ion mode. If FALSE, no additional adduct types will be searched.
#' @param processing.algorithm A single character. "Default", "compMS2Miner" or "RMassBank". 
#' @param params.search Parameters for searching and collecting ions from chromatogram files in a list. These parameters define the tolerance window when input metadata is searched. The list must contain following elements:
#' \itemize{
#'  \item{mz_search:}{ Numeric. Absolute mass tolerance in Da.}
#'  \item{ppm_search:}{ Numeric. Absolute mass tolerance in ppm.} 
#'  \item{rt_search":}{ Numeric. Absolute retention time tolerance in second.}
#'  \item{rt_gap:}{ Numeric. Retention time gap in second - when two scans both match with an input structure, they are both recorded as isomeric features of the same identifier if they are separated by a certain retention time gap. Please set it to 10000 if no isomeric feature is picked. This parameter is not used for RMassBank.}
#' }
#' @param params.ms.preprocessing Paremeters for post-processing scans found in chromatogram files in a list. It must contain:
#' \itemize{
#'  \item{normalized:}{ Logical. TRUE if the intensities of extracted spectra need to normalized so that the intensity of highest peak will be 100.}
#'  \item{baseline:}{ Numeric. Absolute intensity threshold that is considered as a mass peak and written into the library.}
#'  \item{relative:}{ Numeric between 0 and 100. Relative intensity threshold of the highest peak in each spectrum, peaks above both absolute and relative thresholds are saved in the library.}
#'  \item{max_peaks:}{ Integer higher than 3. Maximum number of peaks kept per spectrum from the highest peak.}
#'  \item{recalibration:}{ NUmeric. Parameter used by RMassBank. 0 if output is experimental spectra. 1 if output is experimental mass along with annotated formula. 2 if output is the theoritical masses calculated from elemenetal formula.}
#' } 
#' @param params.user A single character. User name who processed the data.
#' 
#' @return
#' \itemize{
#'   \item{complete:}{ Entire spectra library (historical + newly added records) is a list object of two elements: "library$sp" ~ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity; "library$metadata" ~ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on scans detected in the chromatogram files. Following metadata columns are updated/added: FILENAME (which raw data file the scan is isolated), MSLEVEL (1 or 2), TIC, PEPMASS_DEV (ppm error for detected precursor mass) and SCANNUMBER (scan number in raw chromatogram). The last three columns were PARAM_SUBMIT_USER (user name), PARAM_ALGORITHM (algorithm of processing), PARAM_CREATION_TIME (date and time when the MS record was added) and SCANS (unique identifier for each record)}
#'   \item{current:}{ Temporary spectra library that only contains newly added scans.}
#' }
#'
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
#'
#' @examples
#'
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils read.csv
#' @importFrom plyr rbind.fill
#' 
#' @export
#'
#' @examples
#'
#' library(RMassBankData)
#'
#' input_library = NULL # There's no historical spectral library. We create a brand new spectral library here,
#' raw_data_files <- list.files(system.file("spectra", package="RMassBankData"), ".mzML", full.names = TRUE)
#' metadata_file <- list.files(system.file(package = "MergeION"),".csv", full.names = TRUE)
#' 
#' polarity = "Positive"
#' mslevel= 2 # Only MS2 scans are extracted!
#' add.adduct = F # No additional adducts are searched besides M+H 
#' 
#' params.search = list(mz_search = 0.005, ppm_search = 10, rt_seach = 15, rt_gap = 30)
#' params.ms.preprocessing = list(normalized = T, baseline = 1000, relative =0.01, max_peaks = 200, recalibration = 0)
#' params.user = "DANIEL"
#' 
#' processing.algorithm = "Default"
#' lib = library_generator(input_library, raw_data_files, metadata_file, 
#'                        polarity = "Positive", mslevel, add.adduct, processing.algorithm,
#'                        params.search, params.ms.preprocessing, params.user)
#' lib1 = lib$complete
#' save(lib1, file = "test_Default.RData") # Save the library as RData
#' 
#' processing.algorithm = "compMS2Miner"
#' lib = library_generator(input_library, raw_data_files, metadata_file, 
#'                        polarity = "Positive", mslevel, add.adduct, processing.algorithm,
#'                        params.search, params.ms.preprocessing, params.user)
#' lib2 = lib$complete
#' save(lib2, file = "test_compMS2Miner.RData")  # Save the library as RData
#'
#' # Processing with RMassBank:
#' processing.algorithm = "RMassBank"
#' lib = library_generator(input_library, raw_data_files, metadata_file, 
#'                       polarity = "Positive", mslevel, add.adduct, processing.algorithm,
#'                       params.search, params.ms.preprocessing, params.user)
#' lib3 = lib$complete
#' save(lib3, file = "test_RMassBank.RData")  # Save the library as RData
#' 
#' # Processing with RMassBank and only keep the exact fragment masses based on elemental formulas:
#' params.ms.preprocessing$recalibration = 2
#' lib = library_generator(input_library, raw_data_files, metadata_file, 
#'                      polarity = "Positive", mslevel, add.adduct, processing.algorithm,
#'                      params.search, params.ms.preprocessing, params.user)
#" lib4 = lib$complete
#' save(lib4, file = "test_RMassBank_bis.RData")  
#' 
library_generator<-function(input_library = NULL, raw_data_files = NULL, metadata_file = NULL, 
                  polarity = c("Positive", "Negative"), mslevel = c(1, 2), add.adduct = TRUE,
                  processing.algorithm = c("Default", "compMS2Miner", "RMassBank"),
                  params.search = list(mz_search = 0.005, ppm_search = 10, rt_seach = 15, rt_gap = 30), 
                  params.ms.preprocessing = list(normalized = T, baseline = 1000, relative =0.01, max_peaks = 200, recalibration = 0),
                  params.user = ""){

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  #####################################
  ### Check general function inputs:###
  #####################################

  if (is.character(input_library)){
    if (input_library!=""){
      if (file_ext(input_library)!="mgf" & file_ext(input_library)!="RData" & file_ext(input_library)!="msp"){
        stop("The input library must be mgf, msp or RData format!")
   }}}
  
  if (!is.vector(raw_data_files)){
    stop("Please provide a list of chromatogram files!")}

  if (!all(file_ext(raw_data_files) %in% c("mzML", "mzXML", "cdf", "CDF"))){
    stop("Chromatogram files must be in mzML, mzXML or cdf format!")}

  if (!is.character(metadata_file)){
    stop("You must provide the name of the csv or excel file that contains targeted metabolic features!")
  }
  
  if (length(metadata_file)>1){
      stop("Metadata must be written in one csv, txt or dat file!")
  } 
  
  if (file_ext(metadata_file)!="csv" & file_ext(metadata_file)!="txt" & file_ext(metadata_file)!="dat"){
    stop("Metadata must be written in one csv, txt or dat file!")
  } 
  
  if (length(polarity)!=1){
    stop("Polarity must be either Positive or Negative")
  }
  
  if (!(polarity %in% c("Positive", "Negative"))){
    stop("Polarity must be either Positive or Negative")
  }
  
  if (!(2 %in% mslevel)){
    stop("2 must in mslevel!")
  }
  
  if (length(processing.algorithm)!=1){
    stop("Please choose one processing algorithm!")
  }
  
  if (!(processing.algorithm %in% c("Default", "compMS2Miner", "RMassBank"))){
    processing.algorithm = "Default"
    message("Default algorithm is set, otherwise please set processing.algorithm as compMS2Miner or RMassBank")
  }
  
  if (!(processing.algorithm == "RMassBank") & params.ms.preprocessing$recalibration!=0){
    message("Recalibration or Elemental formula calculation can be performed with RMassBank algorithm only!")
  }
  
  ###############################
  ### Read and check old library:
  ###############################
  
  if (is.character(input_library)){
    if (file_ext(input_library)=="mgf"){
      old_lib = readMGF2(input_library)}
    if (file_ext(input_library)=="RData"){
      old_lib = load_object(input_library)
    if (file_ext(input_library)=="msp"){
      old_lib = readMSP2(input_library)}
    }
  } else {old_lib=input_library}
  
  if (!is.null(old_lib)){
    if (length(old_lib)==2 & "complete" %in% names(old_lib)){
      old_lib = old_lib$complete}
    spectrum_list=old_lib$sp
    metadata=old_lib$metadata
    max_scan = max(as.numeric(metadata$SCANS)) # Maximal scan
    NN = length(spectrum_list)
  } else {
    spectrum_list = list()
    metadata = c()
    max_scan = 0
    NN= 0 
  }
  
  # Old library size:
  NN0 = NN

  #############################
  ### Load and edit metadata###
  #############################
  
  ref = read.csv(metadata_file,sep=";",dec=".",header=T)
  if (ncol(ref)==1){ref = read.csv(metadata_file,sep=",",dec=".",header=T)}  
  if (ncol(ref)==1){ref = read.csv(metadata_file,sep="\t",dec=".",header=T)}  
  if (ncol(ref)<2){stop("Input metadata format not valid!")}

  if (!("ID" %in% colnames(ref))){stop("Input metadata must contain a column ID!")}
  
  if (processing.algorithm=="Default"){
    if (!("PEPMASS" %in% colnames(ref))){stop("To use Default algorithm, PEPMASS must be provided! You could set it to N/A if SMILES is known!")}
  }
  
  if (processing.algorithm=="compMS2Miner"){
    if (!("PEPMASS" %in% colnames(ref))){stop("To use compMS2Miner algorithm, PEPMASS must be provided!")}
    if (!("RT" %in% colnames(ref))){stop("To use compMS2Miner algorithm, retention time (RT) must be provided!")}
  }
  
  if (processing.algorithm=="RMassBank"){
    if (!("SMILES" %in% colnames(ref))){stop("To use RMassBank algorithm, SMILES must be provided!")}
    if (!("FILENAME" %in% colnames(ref))){stop("To use RMassBank algorithm, please provide LC-MS file name corresponding to each compound!")}
  }
  
  target.ref = metadata_editor(ref, processing.algorithm, polarity, add.adduct)
  
  if (is.null(target.ref)){
    stop("No valid metadata available!")
  }
  if (nrow(target.ref)==0){
    stop("No valid metadata available!")
  }

  #################
  ### MergeION#####
  #################

  if (processing.algorithm=="Default"){
  
    FF = length(raw_data_files)
    temp_metadata = c() # Temporary metadata
  
    for (ff in 1:FF){
      
      print(ff)
      # Extract MS2 scans:
      
      if (2 %in% mslevel){
    
        dat2 = process_MS2(raw_data_files[ff], ref = target.ref, 
                           rt_search = params.search$rt_seach, rt_gap = params.search$rt_gap, ppm_search = params.search$ppm_search, mz_search = params.search$mz_search, 
                           baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, normalized = params.ms.preprocessing$normalized)
        
        LL2 = length(dat2$sp) # Added library size
        
        if (LL2>0){
            temp.ref = dat2$ref_MS2 # Filter again metadata data for MS1 searcch
            new_scans2 = (max_scan+1):(max_scan+LL2)
            max_scan = max_scan+LL2
            metadata2 = cbind.data.frame(dat2$metadata, 
                                         PARAM_SUBMIT_USER = params.user, 
                                         PARAM_ALGORITHM = "MergeION",
                                         PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans2)
            for (n in 1:LL2){spectrum_list[[NN+n]]=dat2$sp[[n]]} # Update spectrum list
            temp_metadata = rbind(temp_metadata, metadata2)
            NN=NN+LL2
    }}

     # Extract MS1 scans:

     if (1 %in% mslevel){ # We search MS1 only for compounds that are fragmented to provide isotopic pattern knowledge
        
        dat1 = process_MS1(raw_data_files[ff], ref = temp.ref,
                           rt_search = params.search$rt_seach, rt_gap = params.search$rt_gap, ppm_search = params.search$ppm_search, mz_search = params.search$mz_search, 
                           baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, normalized = params.ms.preprocessing$normalized)
        
        LL1= length(dat1$sp) # Added library size
        if (LL1>0){
          new_scans1 = (max_scan+1):(max_scan+LL1)
          max_scan = max_scan+LL1
          metadata1 = cbind.data.frame(dat1$metadata, 
                                       PARAM_SUBMIT_USER = params.user, 
                                       PARAM_ALGORITHM = "MergeION",
                                       PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans1)
          for (n in 1:LL1){spectrum_list[[NN+n]]=dat1$sp[[n]]} # Update spectrum list
          temp_metadata =  rbind.data.frame(temp_metadata, metadata1)
          NN=NN+LL1
        }}
     }
  
    metadata = rbind.fill(metadata, temp_metadata)
  }
  
  #################
  ### compMS2######
  #################
  
  if (processing.algorithm=="compMS2Miner"){
    
    FF = length(raw_data_files)
    temp_metadata = c() # Temporary metadata
    
    include.MS1 = FALSE
    if (1 %in% mslevel){include.MS1 = TRUE}
    
    for (ff in 1:FF){
      
      print(ff)
      
      dat12 = process_compMS2Miner(raw_data_files[ff], ref = target.ref, polarity = polarity, include.MS1 = include.MS1,
         rt_search = params.search$rt_seach, rt_gap = params.search$rt_gap, ppm_search = params.search$ppm_search, mz_search = params.search$mz_search, 
         baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, normalized = params.ms.preprocessing$normalized)
    
      LL12 = length(dat12$sp) # Added library size
      
      if (LL12>0){
        new_scans12 = (max_scan+1):(max_scan+LL12)
        max_scan = max_scan+LL12
        metadata12 = cbind.data.frame(dat12$metadata, 
                                      PARAM_SUBMIT_USER = params.user,
                                      PARAM_ALGORITHM = "compMS2Miner",
                                      PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans12)
        for (n in 1:LL12){spectrum_list[[NN+n]]=dat12$sp[[n]]} # Update spectrum list
        temp_metadata = rbind(temp_metadata, metadata12)
        NN=NN+LL12
      }}
    metadata = rbind.fill(metadata, temp_metadata)
  }
  
  ###################
  ### RMassBank######
  ###################
  
  if (processing.algorithm=="RMassBank"){
    
    if (!("SMILES" %in% colnames(ref))){stop("To use RMassBank algorithm, valid SMILES code must be provided!")}
        
    FF = length(raw_data_files)
    temp_metadata = c() # Temporary metadata
    
    include.MS1 = FALSE
    if (1 %in% mslevel){include.MS1 = TRUE}

    for (ff in 1:FF){
      
      dat12 = process_RMassBank(raw_data_files[ff], ref = target.ref, polarity = polarity, include.MS1 = include.MS1,
                                rt_search = params.search$rt_seach, ppm_search = params.search$ppm_search, 
                                baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, 
                                recalibration = params.ms.preprocessing$recalibration, normalized = params.ms.preprocessing$normalized)
      
      LL12 = length(dat12$sp) # Added library size
      
      if (LL12>0){
        new_scans12 = (max_scan+1):(max_scan+LL12)
        max_scan = max_scan+LL12
        metadata12 = cbind.data.frame(dat12$metadata, 
                                      PARAM_SUBMIT_USER = params.user,
                                      PARAM_ALGORITHM = "RMassBank",
                                      PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans12)
        for (n in 1:LL12){spectrum_list[[NN+n]]=dat12$sp[[n]]} # Update spectrum list
        temp_metadata = rbind(temp_metadata, metadata12)
        NN=NN+LL12
    }}
    metadata = rbind.fill(metadata, temp_metadata)
    
    unlink("mysettings.ini")
    unlink("Compoundlist.csv")
  }
  
  ####################
  ### Return results:
  ####################

  NN = length(spectrum_list)

  library = list()
  library$sp = spectrum_list
  library$metadata = metadata

  library_current = list()
  library_current$sp = spectrum_list[(NN0+1):NN]
  library_current$metadata = metadata[(NN0+1):NN,]

  return(list(complete = library, current = library_current))
}

######################
### Internal function:
######################

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}
