#' Digitizing and networking LC-MS/MS data
#'
#' The function proposes three data processing algorithms to pick up MS1/MS2 scans from DDA or targeted mode LC-MS/MS data, merge them into a spectral library and create a spectral similarity-based molecular network.
#' 
#' @param input_library Character or a list object. If character, name of the existing library into which new scans are added, the file extension must be mgf, msp or RData; please set to NULL if the new library has no dependency with previous ones.
#' @param lcms_files A character vector of LC-MS/MS file names from which scans are extracted. All files must have be in centroid-mode with mzML, mzXML or cdf extension!
#' @param metadata_file A single character, NULL object or data frame. If it is character, it should be the metadata file name. The file should be tab, comma or semi-colon separated txt, dat or csv format. For all algorithms, the metadata must contain the column "ID" - a unique structure identifier. The column PEPMASS (targeted precursor mass) must be provided for Default and compMS2Miner. The column RT (targeted retention time in min) must be provided for compMS2Miner and optional for MergeION and RMassBank. Please include the column SMILES (structure identifier) for RMassBank algorithm. If RMassBank is used, the column FILENAME (chromatogram file with mzML, mzXML or cdf extension) must be provided for each compound telling the algorithm from which file compound can be found. Column FILENAME is optional for Default and compMS2Miner. Column ADDUCT is optional for all algorithms, if not provided, all input will be considered as M+H or M-H depending on polarity. Please specify the adduct type if metadata contains both positive and negative ions. If metadata is NULL and lcms files are acquired in DDA mode, an automated feature screening is performed for fragmented masses. Masses and retention times of these features are used for spectral library generation and molecular networking. 
#' @param polarity A single character. Either "Positive" or "Negative". Ion mode of LC-MS/MS files. 
#' @param mslevel A numeric vector. 1 or 2 or c(1,2). 2 if MS2 scans are extracted, 1 if isotopic pattern of the precursor mass in the MS1 scan is extracted. c(1,2) if both MS1 and MS2 scans are extracted. Note: High-quality isotopic patterns in MS1 scans are useful for determining precursor formula!
#' @param add.adduct Logical. If TRUE, additional adduct types will be calculated based on precursor masses of "M+H" and "M-H" adducts in the input metadata: "M+2H", "M+Na","M+K","M+NH4","M+" will be searched for positive ion mode, "M+COO-", "M+Cl" and "M+CH3COO-" for negative ion mode. If FALSE, no additional adduct types will be searched.
#' @param adductType. User-specified adduct type, default is NULL. Set `add.adduct` to TRUE and specify `adductType` to fiter records limited to `adductType` before appending the additional adduct types.
#' @param processing.algorithm A single character. "Default", "compMS2Miner" or "RMassBank". 
#' @param params.search Parameters for searching and collecting ions from chromatogram files in a list. These parameters define the tolerance window when input metadata is searched. The list must contain following elements:
#' \itemize{
#'  \item{mz_search:}{ Numeric. Absolute mass tolerance in Da.}
#'  \item{ppm_search:}{ Numeric. Absolute mass tolerance in ppm.} 
#'  \item{rt_search":}{ Numeric. Absolute retention time tolerance in second.}
#'  \item{rt_gap:}{ Numeric. Retention time gap in second - when two scans both match with an input structure, they are both recorded as isomeric features of the same identifier if they are separated by a certain retention time gap. Please set it to 10000 if no isomeric feature is picked. This parameter is not used for RMassBank.}
#' }
#' @param params.ms.preprocessing Parameters for pre-processing scans found in chromatogram files in a list. It must contain:
#' \itemize{
#'  \item{normalized:}{ Logical. TRUE if the intensities of extracted spectra need to normalized so that the intensity of highest peak will be 100.}
#'  \item{baseline:}{ Numeric. Absolute intensity threshold that is considered as a mass peak and written into the library.}
#'  \item{relative:}{ Numeric between 0 and 100. Relative intensity threshold of the highest peak in each spectrum, peaks above both absolute and relative thresholds are saved in the library.}
#'  \item{max_peaks:}{ Integer higher than 3. Maximum number of peaks kept per spectrum from the highest peak.}
#'  \item{recalibration:}{ NUmeric. Parameter used by RMassBank. 0 if output is experimental spectra. 1 if output is experimental mass along with annotated formula. 2 if output is the theoritical masses calculated from elemenetal formula.}
#' } 
#' @param params.consensus Parameters for generating consensus scans that combine spectra of the same compound ID
#' \itemize{ 
#'  \item{consensus:}{ Logical. TRUE if consensus spectra are generated }
#'  \item{consensus_method:}{ Character. Method for merging library "duplicates" by compound IDs.}
#'    \enumerate{
#'     \item{consensus:}{ Default method for generating generated by merging spectra of the same compound ID. All peaks were kept, similar fragments were aligned by averaging m/z and intensity.}
#'     \item{common_peaks:}{ Peaks detected in ALL duplicated spectra were kept and averaged.}
#'     \item{most_recent:}{ The most recent record was kept if duplicates are detected.}}
#'  \item{consensus_window}{ m/z window (in Dalton) for spectra alignment, only used when method = "consensus" or "common_peaks". To generate consensus spectra, mass peaks in different spectra within the mass window were aligned by averaging their mass values and intensities. The metadata was kept only for the most recent spectrum.}
#' }
#' @param params.network Parameters for networking consensus spectra library into a molecular network
#' \itemize{
#'   \item{network:}{ Logical. TRUE if a network is built for consensus spectral library}
#'   \item{similarity_method:}{Characeter.Similarity metrics for networking and spectral library search. Must be "Matches", "Dot", "Cosine", "Spearman", "MassBank" or "NIST".}
#'   \item{min_frag_match:}{ Integer. Minimum number of common fragment ions (or neutral losses) that are shared to be considered for spectral similarity evaluation. We suggest setting this value to at least 6 for statistical meaningfulness.}
#'   \item{min_score:}{ Numeric between 0 and 1. Minimum similarity score to connect two nodes in the library annotate an unknown feature with spectral library or to connect two unknown features because they are similar. It does NOT affect method = "Matches".}
#'   \item{topK:}{ Integer higher than 0. For networking, the edge between two nodes are kept only if both nodes are within each other's TopK most similar nodes. For example, if this value is set at 20, then a single node may be connected to up to 20 other nodes. Keeping this value low makes very large networks (many nodes) much easier to visualize. We suggest keeping this value at 10.}
#'   \item{max_comp_size:}{ Numeric between 0 and 200. Maximum size of nodes allowed in each network component. Default value is 100. Network component = Cluster of connected node. Set to 0 if no limitation on componet size.}
#'   \item{reaction_type:}{ Character. Either "Metabolic" and "Chemical". Type of transformation list used to annotate mass difference between connected features in molecular network.}
#'   \item{use_reaction:}{ Boolean. TRUE if keep only edges whose mass difference can be annotated to known metabolic or chemical reactions.}
#' }
#' @param params.user A list of additional parameters.
#' \itemize{
#'   \item{sample_type:}{ Character. Type of LCMS samples added to the spectral library e.g. plasma, standards...}
#'   \item{user_name:}{ Character. User name who process the batch of lc-ms files.}
#'   \item{comments:}{ Character. Additional comments about the samples added.}
#' }
#' 
#' @return
#' \itemize{
#'   \item{complete:}{ Entire spectra library (historical + newly added records) is a list object of two elements: "library$sp" ~ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity; "library$metadata" ~ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on scans detected in the chromatogram files. Following metadata columns are updated/added: FILENAME (which raw data file the scan is isolated), MSLEVEL (1 or 2), TIC, PEPMASS_DEV (ppm error for detected precursor mass) and SCANNUMBER (scan number in raw chromatogram). The last three columns were PARAM_ALGORITHM (algorithm of processing), PARAM_CREATION_TIME (date and time when the MS record was added) and SCANS (unique identifier for each record)}
#'   \item{consensus:}{ Consensus spectral library by merging MS/MS spectra with the same ID.}
#'   \item{network:}{ Consensus spectral library transformed into a molecular network based on MS/MS spectral similarity.}
#' }
#'
#' @author Youzhong Liu, \email{YLiu186@ITS.JNJ.com}
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
#' lcms_files <- list.files(system.file("spectra", package="RMassBankData"), ".mzML", full.names = TRUE)
#' metadata_file <- list.files(system.file(package = "MergeION"),".csv", full.names = TRUE)
#' 
#' polarity = "Positive"
#' mslevel= 2 # Only MS2 scans are extracted!
#' add.adduct = F # No additional adducts are searched besides M+H 
#' 
#' params.search = list(mz_search = 0.005, ppm_search = 10, rt_search = 15, rt_gap = 30)
#' params.ms.preprocessing = list(normalized = T, baseline = 1000, relative =0.01, max_peaks = 200, recalibration = 0)
#' 
#' # Building a spectral library with default (SmartION) algorithm by simply gathering scans that matched with metadata:
#' params.user = list(sample_type = "RMassBank data", user_name = "daniel", comments = "default algorithm, without building a consensus library")
#' processing.algorithm = "Default"
#' lib = library_generator(input_library, lcms_files, metadata_file, 
#'                         polarity = "Positive", mslevel, add.adduct, processing.algorithm,
#'                         params.search, params.ms.preprocessing, params.user = params.user)
#' lib1 = lib$complete
#' save(lib1, file = "test_default_complete.RData") # Save the library as RData
#' 
#' # Building a spectral library with compMS2Miner algorithm and generating a consensus spectral library
#' processing.algorithm = "compMS2Miner"
#' params.consensus = list(consensus = T, consensus_method = "consensus", consensus_window = 0.02)
#' params.user = list(sample_type = "RMassBank data", user_name = "daniel", comments = "compMS2Miner algorithm, building a consensus library")
#' lib = library_generator(input_library, lcms_files, metadata_file, 
#'                        polarity = "Positive", mslevel, add.adduct, processing.algorithm,
#'                        params.search, params.ms.preprocessing, params.consensus, params.user = params.user)
#' lib2 = lib$consensus
#' save(lib2, file = "test_compMS2Miner_consensus.RData")  # Save the library as RData
#'
#' # Building a spectral library with RMassBank algorithm (recalibration based on elemental formula annotation), creating consensus spectral library and building a molecular network based on the consensus library
#' processing.algorithm = "RMassBank"
#' params.ms.preprocessing = list(normalized = T, baseline = 1000, relative =0.01, max_peaks = 200, recalibration = 2)
#' params.network = list(network = T, similarity_method = "Cosine", min_frag_match = 6, min_score = 0.6, max_comp_size = 100, topK = 10, reaction_type = "Chemical", use_reaction = F)
#' lib3 = library_generator(input_library, lcms_files, metadata_file, 
#'                       polarity = "Positive", mslevel, add.adduct, processing.algorithm,
#'                       params.search, params.ms.preprocessing, params.consensus, params.network, params.user = params.user)
#' save(lib3, file = "test_RMassBank_consensus_network.RData")  # Save the library as RData
#' 
library_generator<-function(input_library = NULL, lcms_files = NULL, metadata_file = NULL, 
                  polarity = c("Positive", "Negative")[1], mslevel = c(1, 2), add.adduct = TRUE, adductType = NULL,
                  processing.algorithm = c("Default", "compMS2Miner", "RMassBank")[1],
                  params.search = list(mz_search = 0.01, ppm_search = 10, rt_search = 15, rt_gap = 30), 
                  params.ms.preprocessing = list(normalized = TRUE, baseline = 1000, relative = 0.1, max_peaks = 200, recalibration = 0),
                  params.consensus = list(consensus = FALSE, consensus_method = c("consensus", "consensus2", "common_peaks","most_recent")[1], consensus_window = 0.02),
                  params.network = list(network = FALSE, similarity_method = "Cosine", min_frag_match = 6, min_score = 0.6, topK = 10, max_comp_size = 100, reaction_type = "Metabolic", use_reaction = FALSE),
                  params.user = list(sample_type = "", user_name = "", comments = "")){

  options(stringsAsFactors = FALSE)
  options(warn=-1)
  gc()
  
  output_library = NULL
  
  #####################################
  ### Check general function inputs:###
  #####################################

  if (is.null(input_library) & is.null(lcms_files)){
    stop("Please at least provide input spectral library or a vector of lc ms files!")
  }
  
  old_lib = old_consensus = old_network = NULL
  
  if (!is.null(input_library)){
    old_lib = library_reader(input_library)$complete
    old_consensus = input_library$consensus
    old_network  = input_library$network
  }
  
  if (!is.null(lcms_files)){
    if (!is.vector(lcms_files)){
      stop("Please provide a list of chromatogram files!")
    }
  
    if (!all(toupper(file_ext(lcms_files)) %in% c("MZML", "MZXML", "CDF"))){
      stop("Chromatogram files must be in mzML, mzXML or cdf format!")
    }
  }

  if (!is.null(metadata_file)){
    if (is.list(metadata_file)){
       ref = metadata_file
       metadata_file = "durlach.csv"
    }
    if (!is.character(metadata_file)){
      stop("You must provide the name of the csv or excel file that contains targeted metabolic features!")
    }
    if (length(metadata_file)>1){
      stop("Metadata must be written in one csv, txt or dat file!")
    } 
    if (file_ext(metadata_file)!="csv" & file_ext(metadata_file)!="txt" & file_ext(metadata_file)!="dat"){
      stop("Metadata must be written in one csv, txt or dat file!")
    }
  }
  
  if (is.null(metadata_file)){add.adduct = FALSE}
  
  if (length(polarity)!=1){
    stop("Polarity must be either Positive or Negative")
  }
  
  if (!(polarity %in% c("Positive", "Negative"))){
    stop("Polarity must be either Positive or Negative")
  }
  
  if (length(processing.algorithm)!=1){
    stop("Please choose one processing algorithm!")
  }
  
  if (!(processing.algorithm %in% c("Default", "compMS2Miner", "RMassBank"))){
    processing.algorithm = "Default"
    message("Default algorithm (SMartION) is set, otherwise please set processing.algorithm as compMS2Miner or RMassBank!")
  }
  
  if (processing.algorithm=="RMassBank" & is.null(metadata_file)){
    processing.algorithm = "Default"
    message("Default algorithm (SMartION) is set because RMassBank only works if struture metadata is provided!")
  }

  if (!(processing.algorithm == "RMassBank") & params.ms.preprocessing$recalibration!=0){
    message("Recalibration or Elemental formula calculation can be performed with RMassBank algorithm only!")
  }
  
  if (length(mslevel)==1){
    if (mslevel==1){
      processing.algorithm=="Default"
      message("Processing algorithm is set to default since only MS1 scan is extracted!")
    }
  }
  
  if (!params.consensus$consensus & params.network$network){
    params.consensus$consensus = TRUE
  }
  
  if (!params.consensus$consensus){
    message("Warning! Spectral library search impossible without consensus spectral library generation!")
  }
  
  #######################
  ### Check old library:#
  #######################

  if (!is.null(old_lib)){
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
  
  FF = length(lcms_files)
  target.ref = NULL
  
  if (is.null(metadata_file) & FF>0){
    ref = c()
    for (x in 1:length(lcms_files)){
      tmp_ref = process_dda(lcms_files[x], polarity= polarity,  ppm_search = params.search$ppm_search, rt_search = params.search$rt_search, baseline = params.ms.preprocessing$baseline)
      ref = rbind.data.frame(ref, tmp_ref)
    }
    target.ref = ref
    if (nrow(target.ref)==0){target.ref= NULL}
  }
  
  if (!is.null(metadata_file)){
    
    if (metadata_file=="durlach.csv"){
    } else {
      ref = read.csv(metadata_file,sep=";",dec=".",header=T)
      if (ncol(ref)==1){ref = read.csv(metadata_file,sep=",",dec=".",header=T)}  
      if (ncol(ref)==1){ref = read.csv(metadata_file,sep="\t",dec=".",header=T)}  
    }
      
  if (ncol(ref)<2){stop("Input metadata format not valid!")}
    
  if (!("ID" %in% colnames(ref))){stop("Input metadata must contain a column ID!")}
  
  if ("FILENAME" %in% colnames(ref)){
      valid = which(basename(ref$FILENAME) %in% basename(lcms_files))
      ref = ref[valid,,drop = FALSE]
      valid = which(basename(lcms_files) %in% ref$FILENAME)
      lcms_files = lcms_files[valid]
      
      if (length(lcms_files)==0){stop("FILENAME in metadata does not match with raw lcms file name!")}
  }
    
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
    target.ref = process_metadata(ref, processing.algorithm, polarity, add.adduct, adductType)
    if (nrow(target.ref)==0){target.ref= NULL}
  }
  
  if (is.null(target.ref)){FF=0}
  
  if (FF==0){} else{
  
  ################
  ### SmartION####
  ################
  
  if (processing.algorithm=="Default"){
  
    temp_metadata = c() # Temporary metadata
  
    for (ff in 1:FF){
      
      print(ff)
      temp.ref = target.ref
      
      # Extract MS2 scans:
      
      if (2 %in% mslevel){
        
        dat2 = process_SmartMS2(lcms_files[ff], ref = target.ref, 
            rt_search = params.search$rt_search, rt_gap = params.search$rt_gap, ppm_search = params.search$ppm_search, mz_search = params.search$mz_search, 
            baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, normalized = params.ms.preprocessing$normalized)
        
        LL2 = length(dat2$sp) # Added library size
      
        if (LL2>0){
            temp.ref = dat2$ref_MS2 # Filter again metadata data for MS1 searcch
            new_scans2 = (max_scan+1):(max_scan+LL2)
            max_scan = max_scan+LL2
            metadata2 = cbind.data.frame(dat2$metadata, 
                            PARAM_FLAG = 0,
                            PARAM_SUBMIT_USER = params.user$user_name, 
                            PARAM_SAMPLE_TYPE = params.user$sample_type, 
                            PARAM_COMMENTS = params.user$comments, 
                            PARAM_ALGORITHM = "SmartION",
                            PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans2)
            for (n in 1:LL2){spectrum_list[[NN+n]]=dat2$sp[[n]]} # Update spectrum list
            temp_metadata = rbind(temp_metadata, metadata2)
            NN=NN+LL2
    }} 

     # Extract MS1 scans:

     if (1 %in% mslevel){ # We search MS1 only for compounds that are fragmented to provide isotopic pattern knowledge

       print(ff)
      
       dat1 = process_SmartMS1(lcms_files[ff], ref = temp.ref,
                  rt_search = params.search$rt_search, rt_gap = params.search$rt_gap, ppm_search = params.search$ppm_search, mz_search = params.search$mz_search, 
                  baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, normalized = params.ms.preprocessing$normalized)

        LL1= length(dat1$sp) # Added library size
        if (LL1>0){
          new_scans1 = (max_scan+1):(max_scan+LL1)
          max_scan = max_scan+LL1
          metadata1 = cbind.data.frame(dat1$metadata, 
                            PARAM_FLAG = 0,
                            PARAM_SUBMIT_USER = params.user$user_name, 
                            PARAM_SAMPLE_TYPE = params.user$sample_type, 
                            PARAM_COMMENTS = params.user$comments, 
                            PARAM_ALGORITHM = "SmartION",
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
    
    temp_metadata = c() # Temporary metadata
    
    include.MS1 = FALSE
    if (1 %in% mslevel){include.MS1 = TRUE}
    
    for (ff in 1:FF){
      
      dat12 = process_compMS2Miner(lcms_files[ff], ref = target.ref, polarity = polarity, include.MS1 = include.MS1,
         rt_search = params.search$rt_search, rt_gap = params.search$rt_gap, ppm_search = params.search$ppm_search, mz_search = params.search$mz_search, 
         baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, normalized = params.ms.preprocessing$normalized)
    
      LL12 = length(dat12$sp) # Added library size
      
      if (LL12>0){
        new_scans12 = (max_scan+1):(max_scan+LL12)
        max_scan = max_scan+LL12
        metadata12 = cbind.data.frame(dat12$metadata, 
                                      PARAM_SUBMIT_USER = params.user$user_name, 
                                      PARAM_SAMPLE_TYPE = params.user$sample_type, 
                                      PARAM_COMMENTS = params.user$comments, 
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
        
    FF = length(lcms_files)
    temp_metadata = c() # Temporary metadata
    
    include.MS1 = FALSE
    if (1 %in% mslevel){include.MS1 = TRUE}

    for (ff in 1:FF){
      
      dat12 = process_RMassBank(lcms_files[ff], ref = target.ref, polarity = polarity, include.MS1 = include.MS1,
                            rt_search = params.search$rt_search, ppm_search = params.search$ppm_search, 
                            baseline = params.ms.preprocessing$baseline, relative = params.ms.preprocessing$relative, max_peaks = params.ms.preprocessing$max_peaks, 
                            recalibration = params.ms.preprocessing$recalibration, normalized = params.ms.preprocessing$normalized)
      
      LL12 = length(dat12$sp) # Added library size
      
      if (LL12>0){
        new_scans12 = (max_scan+1):(max_scan+LL12)
        max_scan = max_scan+LL12
        metadata12 = cbind.data.frame(dat12$metadata, 
                                      PARAM_SUBMIT_USER = params.user$user_name, 
                                      PARAM_SAMPLE_TYPE = params.user$sample_type, 
                                      PARAM_COMMENTS = params.user$comments, 
                                      PARAM_ALGORITHM = "RMassBank",
                                      PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans12)
        for (n in 1:LL12){spectrum_list[[NN+n]]=dat12$sp[[n]]} # Update spectrum list
        temp_metadata = rbind(temp_metadata, metadata12)
        NN=NN+LL12
    }}
    metadata = rbind.fill(metadata, temp_metadata)
    
    smile_corrected = sapply(metadata$SMILES, function(x) strsplit(x, "\\.")[[1]][1])
    metadata$SMILES = as.character(smile_corrected)
    
    unlink("mysettings.ini")
    unlink("Compoundlist.csv")
    }
  }
  
  ###################################
  ### Combine and validate results:##
  ###################################

  NN = length(spectrum_list)

  if (NN>0){
    library_complete = list()
    library_complete$sp = spectrum_list
    library_complete$metadata = metadata
    library_complete = remove_blanks(library_complete)
    output_library = library_reader(library_complete)
    NN = nrow(output_library$complete$metadata)
  }
  
  ##############################
  ### Post-processing library###
  ##############################
  
  if (NN>1 & params.consensus$consensus & is.null(old_consensus$metadata)){
    
    library_consensus = process_consensus(library_complete, params.consensus$consensus_method, params.consensus$consensus_window, 
                  params.ms.preprocessing$relative, params.ms.preprocessing$max_peaks)
    
    output_library = library_reader(library_consensus)
  }
  
  if (NN>1 & !is.null(old_consensus$metadata)){
    output_library = list(complete = library_complete, consensus = old_consensus, network = old_network)
  }
  
  NN = nrow(output_library$consensus$metadata)

  if (!is.null(NN)){
    
    if (NN>1 & params.consensus$consensus){
      
    library_network = process_lib2network(output_library, networking = params.network$network, polarity = polarity, 
        params.search = list(mz_search = params.consensus$consensus_window, ppm_search = params.search$ppm_search),
        params.similarity = list(method = params.network$similarity_method, min.frag.match = params.network$min_frag_match, min.score = params.network$min_score),
        params.network = list(topK = params.network$topK, max.comp.size = params.network$max_comp_size, reaction.type = params.network$reaction_type, use.reaction = params.network$use_reaction))
      
    output_library = library_reader(library_network)
  }}
 
  return(output_library)
}

#############################
######Internal function######
#############################

remove_blanks<-function(library1){
  
  # The function remove emptys record from library
  
  sp= library1$sp
  metadata = library1$metadata
  
  valid = which(!sapply(sp, is.null))
  sp = sp[valid]
  metadata = metadata[valid,,drop=FALSE]
  
  valid = which(sapply(sp, nrow)>0)
  sp = sp[valid]
  metadata = metadata[valid,,drop=FALSE]
  
  valid = which(sapply(sp, function(x) x[1,1])>0)
  sp = sp[valid]
  metadata = metadata[valid,,drop=FALSE]
  
  new_library = list(metadata= metadata, sp = sp)
  
  return(new_library)
}
