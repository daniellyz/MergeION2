#' Generating spectral library from raw LC-MS/MS chromatograms
#'
#' The function picks up targeted MS1/MS2 scans and merge them into a spcetral library (new or existing). The raw LC-MS/MS files must be centroid-mode mzML, mzMXL or mzData. Ions are selected based on m/z (and retention time) specified in the metadata (recommended) or by automatic peak picking in XCMS package.
#'
#' @param raw_data_files A character vector of file names of chromatograms from which scans are extracted. All files must have be in centroid-mode with mzML or mzMXL extension!
#' @param metadata_file A single character or NULL. If provided, it must be the metadata file name with csv extension. The first six columns of the metadata must be (in order): "PEPMASS" (precursor masses that we want to find in chromatograms), "RT" (retention time of metabolic features to be found, in minute, please put it to N/A if unknown), "IONMODE" (must be "Positive" or "Negative"),"ADDUCT" (precursor ion adduct type, must be one of "M+H","M+Na","M+K","M-H" and "M+Cl"), "CHARGE" (charge number, please keep it at 1) and "ID" (A unique identifier for targeted compounds in spectral library). If NULL, please set MS1.screener = TRUE, a non-targeted feature screening will be performed using centWave from XCMS. In current release, this functionality only works when all input files are acquired on the same instrument and from the same ion mode, and they must all contain MS1 scans.
#' @param mslevel Must be 1 (if only MS1 scans/isotopic patterns of targeted m/z are extracted), 2 (if only MS2 scans are extracted) or c(1,2) (if both MS1 and MS2 scans are extracted). Note: Isotopic patterns in MS1 scans are useful for determining precursor formula!
#' @param MS1.screener Logical. TRUE if centWave algorithm from XCMS is used to detect LC-MS feature peaks based on MS1 scans in the data. The detected feature peaks are used as metadata to assist the generation of spectral library.
#' @param MS2_type  A single character ("DDA" or "Targeted") if all raw_dat_files are acquired in the same mode; A character vector precising the acquisition mode of each file in raw_data_files (e.g. c("DDA","Targeted","DDA"))
#' @param adduct_type Vector of character. Adduct types of ions considered. Its elements must be among "Default","M+H","M+Na","M+K","M+NH4","M-H" and "M+Cl". No additional ion species will be calculated if "Default".
#' @param max.charge Integer. Maximal charge number. Must be a positive integer e.g. 2 if +1, +2 (or -1, -2) ions are consired.
#' @param isomers Logical. TRUE if isomers are kept (scans with same precursor mass but with difference in retention time higher than 2 * rt_search). If FALSE, only the isomer with highest TIC is kept.
#' @param rt_search Retention time search tolerance (in second) for targeted RT
#' @param ppm_search m/z search tolerance (in ppm) for targeted m/z
#' @param baseline Numeric. The absolute intensity threshold) that is considered as a mass peak and written into the library. Peaks above both absolute and relative thresholds are saved in the library.
#' @param relative Numeric between 0 and 100. The relative intensity threshold of the highest peak in each spectrum). Peaks above both absolute and relative thresholds are saved in the library
#' @param snthrehold Numeric higher than 1. Only used when MS1.screener = TRUE. Parameter used by centWave in XCMS to define chromatogram peaks.
#' @param normalized Logical. TRUE if the intensities of extracted spectra need to normalized so that the intensity of highest peak will be 100
#' @param user Character. Name or ID of the user(s) that created or updated the library.
#' @param write_files Logical. TRUE if user wishes to write the mgf and metadata (txt) file in the folder
#' @param input_library Character or library object. If character, name of the library into which new scans are added, the file extension must be mgf; please set to empty string "" or NULL if the new library has no dependency with previous ones.
#' @param output_library Character.Name of the output library, the file extension must be mgf
#'
#' @return
#' \itemize{
#'   \item{complete:}{ Entire spectra library (historical + newly added records) is a list object of two elements: "library$sp" ~ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity; "library$metadata" ~ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on actually-detected scans. Following metadata columns are added: FILENAME (which raw data file the scan is isolated), MSLEVEL (1 or 2), TIC, PEPMASS_DEV (ppm error for detected precursor mass) and SCANNUMBER (scan number in raw chromatogram). Parameters used for library generation were appended. The last three columns were PARAM_USER (user name), PARAM_CREATION_TIME (date and time when the MS record was added) and SCANS (unique identifier for each record, unchanged) }
#'   \item{current:}{ Temporary spectra library that only contains newly added scans.}
#'   \item{<ouput_library>:}{ A mgf spectral library file (complete spectralibrary) will be written in user's working directory. It contains both spectra and metadata}
#'   \item{<ouput_library.txt>:}{ Metadata will be written as a tab-seperated .txt file in user's working directory. Users can check this file in excel or open office.}
#' }
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#'
#' @examples
#' ### We download four test data sets:
#'
#' url = "https://zenodo.org/record/2581847/files/"
#' original_files = c("NA_170405_MAS006_10.mzML",
#'                   "TESTMIX2_180504_MAS011_06.mzXML",
#'                   "JNJ42165279_171214_MAS006_14.mzXML",
#'                   "GMP_R601592_150925_MAS006_04.mzXML")
#' download.file(paste0(url,original_files[1]),destfile="MIX1.mzML") # Download and rename the files
#' download.file(paste0(url,original_files[2]),destfile="MIX2.mzXML")
#' download.file(paste0(url,original_files[3]),destfile="JNJ.mzXML")
#' download.file(paste0(url,original_files[4]),destfile="GMP.mzXML")
#'
#' ### We create the first library
#' raw_data_files = c("MIX1.mzML","MIX2.mzXML","JNJ.mzXML")
#' metadata_file = "https://raw.githubusercontent.com/daniellyz/MergeION/master/inst/library_metadata.csv"
#'
#' mslevel = c(1,2)  # Both MS1 and MS2 scans are extracted!
#' MS2_type = c("DDA","DDA","Targeted") # Mode of MS/MS experiment for the three files
#' adduct_type = c("Default") # Only looking for default ion types (ion types provided by users in metadata)
#' max.charge = 1 # Only looking for +1 charged ions
#' isomers = FALSE # If isomers are present, only the peak with higher TIC is extracted.
#'
#' rt_search = 12 # Retention time tolerance (s)
#' ppm_search = 10  # Mass tolerance (ppm)
#' baseline = 1000  # Baseline level 1000 is fixed for 3 datasets.
#' relative = 1 # Relative intensitiy level 1% is fixed. All peaks under both baseline and relative level are considered as noise.
#' normalized = TRUE # The intensities of extracted spectra will be normalized to 100 (the highest peak)
#'
#' write_files = FALSE # The library(mgf) and metadata will not be writen in user's folder
#' input_library = "" # A brand new library, there's no previous dependency
#' output_library = "library_V1.mgf" # Name of the library
#' user_name = "Florian" # User name for uploading
#'
#' library1 = library_generator(raw_data_files, metadata_file, mslevel, MS2_type, adduct_type, max.charge, isomers,
#'                             rt_search, ppm_search, baseline, relative, normalized,
#'                             user = user_name, write_files, input_library, output_library)
#'
#' library1 = library1$complete # Important! We extract the library object. "$complete" for extracting the entire library including historical mass spectra. Here since we create a brand-new library, "library1$complete" and "library1$current" are the same.
#'
#' ### Now we process and add a new data GMP.mzXML in the existing library:
#' raw_data_files = "GMP.mzXML"
#' adduct_type = c("M+H", "M+Na") # Two adduct types are now considered
#' MS2_type = "Targeted"
#' isomers = TRUE # We would like now to record all isomers in the library
#'
#' write_files = TRUE # We want to directly write the library mgf + metadata files
#' input_library = library1
#' output_library = "library_V2.mgf"
#' user_name = "Thomas" # Another user adds records into the library
#' library2 = library_generator(raw_data_files, metadata_file, mslevel, MS2_type, isomers, adduct_type, max.charge,
#'                             rt_search, ppm_search, baseline, relative, normalized,
#'                             user = user_name, write_files, input_library, output_library)
#'
#' # In the end, "library_V2.mgf" should appear in the working directory along with its metadata table (txt files)
#'
#' # Now we check in the newly added scans whether the desired precursor mz is in:
#'
#' tmp_library = library2$current
#' query = library_manager(tmp_library, query = c("PEPMASS = 478.096"), ppm_search = 20)
#' library_visualizer(query)
#'
#' @importFrom MSnbase readMSData rtime tic fData readMgfData precursorMz polarity
#' @importFrom tools file_ext
#' @importFrom utils write.table read.csv menu
#' @importFrom plyr rbind.fill
#' @importFrom pracma findpeaks
#'
#' @export
#'
library_generator<-function(raw_data_files, metadata_file, mslevel = c(1,2), MS1.screener = FALSE,
                            MS2_type = "DDA", isomers = TRUE, adduct_type = "M+H", max.charge = 1,
                            rt_search = 12,ppm_search = 20,
                            baseline = 1000, relative =5, snthreshold = 30, normalized = T,
                            user = "", write_files = TRUE, input_library = "", output_library = ""){

  options(stringsAsFactors = FALSE)
  options(warn=-1)
  FF = length(raw_data_files)
  spectrum_list = list()
  metadata = c()
  NN = 0
  max_scan = 0 #  Highest scan ID
  old_lib = NULL # Previous library
  LL1 = LL2 = 0
  #unlink(output_library)

  ##############################
  ### Check function inputs:
  ##############################

  if (missing(raw_data_files) || (!is.vector(raw_data_files))){
    stop("Please provide a list of chromatogram files!")}

  if (!all(file_ext(raw_data_files) %in% c("mzML","mzXML","mzData"))){
    stop("Chromatogram files must be in mzML, mzXML or mzData format!")}

  if (!is.null(metadata_file)){
    if ((!is.character(metadata_file)) || (length(metadata_file)!=1) || file_ext(metadata_file)!="csv"){
      stop("Metadata must be written in one single csv!")
    } else {
    ref<-read.csv(metadata_file,sep=";",dec=".",header=T,stringsAsFactors = F)}
  } else {ref = NULL}

  if (!all(mslevel %in% c(1,2))){
    stop("mslevel must be 1 or 2!")}

  if (is.null(metadata_file) && !MS1.screener){
    stop("No metadata is available! Please provide metadatafile or set MS1.screener = TRUE!")
  }

  if (length(MS2_type)==1){
    MS2_type = rep(MS2_type,FF)
  } else {
    if (length(MS2_type)!=FF){
      stop("The length of MS2_type must be the same as raw_data_files!")}}

  if (!all(adduct_type %in% c("Default","M+H","M+Na","M+K","M+NH4","M-H","M+Cl"))){
    stop("adduct_type not valid!")}

  if ((max.charge%%1 != 0) || (max.charge<=0)){
    stop("max.charge must be a positive integer!")
  }

  if (length(baseline)==1){
    baseline = rep(baseline,FF)
  } else {
    if (length(baseline)!=FF){
      stop("The length of baseline must be the same as raw_data_files!")}}

  if (length(relative)==1){
    relative = rep(relative,FF)
  } else {
    if (length(relative)!=FF){
      stop("The length of relative must be the same as raw_data_files!")}}

  if (length(snthreshold)==1){
    snthreshold = rep(snthreshold,FF)
  } else {
    if (length(snthreshold)!=FF){
      stop("The length of snthreshold must be the same as raw_data_files!")}}

  MS2_type = match.arg(MS2_type,choices=c("DDA","Targeted"),several.ok = TRUE)

  if (is.character(input_library)){
    if (input_library!=""){
      if (file_ext(input_library)!="mgf"){
        stop("The input library must be mgf format!")
  }}}

  #######################################
  ### Read from metadata and old library:
  #######################################

  if (!is.null(ref)){

    if (is.null(nrow(ref))){ # Only one row...
      labels=names(ref)
      ref=data.frame(matrix(ref,nrow=1))
      colnames(ref)=labels
    }

    ref[,6]=as.character(ref[,6]) # Make sure IDs are characters

    colnames(ref)[1]="PEPMASS"  # Make sure column name is correct!
    colnames(ref)[2]="RT"
    colnames(ref)[3]="IONMODE"
    colnames(ref)[4]="ADDUCT"
    colnames(ref)[5]="CHARGE"
    colnames(ref)[6]="ID"

    for (j in 1:ncol(ref)){  # Change column names to avoid duplicates
      if (colnames(ref)[j] %in% c("FILENAME","MSLEVEL","TIC","PEPMASS_DEV","PEPMASS_DEV","SCAN_NUMBER","SCANS")){
      colnames(ref)[j]=paste0(colnames(ref)[j],"_000")
    }}

    if (ncol(ref)<6){
      stop("Metadata must contain at least 6 columns!")}

    if (!is.numeric(ref$PEPMASS)){
      stop("Precursor masses (PEPMASS) must be numeric!")}

    if (!all(ref$IONMODE %in% c("Positive","Negative"))){
      stop("Ion mode must be Postive or Negative!")}

    if (!all(ref$ADDUCT %in% c("M+H","M+Na","M+K","M-H","M+Cl"))){
      stop("Metadata contain non valid adduct types!")}

    # Calculate multicharge:
    ref = metadata_adduct_editor(ref, adducts = adduct_type, max.charge = max.charge)
  }

  if (is.character(input_library)){
   if (input_library!=""){
      old_lib=readMGF2(input_library)
   }} else {
   if (!is.null(input_library)){
      old_lib=input_library}}

  if (!is.null(old_lib)){
    if (length(old_lib)==2 & "complete" %in% names(old_lib)){
      old_lib = old_lib$complete
    }
    spectrum_list=old_lib$sp
    metadata=old_lib$metadata
    max_scan = max(as.numeric(metadata$SCANS))
    #metadata_items=colnames(metadata)[1:(ncol(metadata)-13)]
    #metadata=metadata[,1:(ncol(metadata)-1)] # Remove LAST COLUMN SCANS!
    NN = length(spectrum_list)
  }

  # Old library size:
  NN0 = NN

  ##################
  ### Batch process:
  ##################

  for (ff in 1:FF){

    # Define metadata before processing each individual file:

    if (MS1.screener){ # Redefine metadata if MS1 screener is on
      targeted.ref = metadata_MS1_screener(raw_data_files[ff], ref = ref,
                          max.charge = max.charge, ppm_search = ppm_search, rt_search = rt_search,
                          baseline = baseline[ff], snthreshold = snthreshold[ff])
    } else {targeted.ref = ref}

    # Extract MS2 scans:

    if (!is.null(targeted.ref)){

      if (2 %in% mslevel){
        dat2 = process_MS2(raw_data_files[ff], targeted.ref, rt_search, ppm_search, MS2_type[ff], isomers, baseline[ff], relative[ff])
        LL2 = length(dat2$sp) # Added library size
        if (LL2>0){
          targeted.ref = dat2$ref_MS2 # Filter again metadata data for MS1 searcch
          new_scans2 = (max_scan+1):(max_scan+LL2)
          max_scan = max_scan+LL2
          metadata2 = cbind.data.frame(dat2$metadata, PARAM_SUBMIT_USER = user, PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans2)
          for (n in 1:LL2){spectrum_list[[NN+n]]=dat2$sp[[n]]} # Update spectrum list
            metadata=rbind.fill(metadata,metadata2) # Update metadata
            NN=NN+LL2
        }}

    # Extract MS1 scans:

      if (1 %in% mslevel){ # We search MS1 only for compounds that are fragmented to provide isotopic pattern knowledge
        dat1 = process_MS1(raw_data_files[ff], targeted.ref, rt_search, ppm_search, isomers, baseline[ff], relative[ff])
        LL1= length(dat1$sp) # Added library size
        if (LL1>0){
          new_scans1 = (max_scan+1):(max_scan+LL1)
          max_scan = max_scan+LL1
          metadata1 = cbind.data.frame(dat1$metadata, PARAM_SUBMIT_USER = user, PARAM_CREATION_TIME = Sys.time(), SCANS = new_scans1)
          for (n in 1:LL1){spectrum_list[[NN+n]]=dat1$sp[[n]]} # Update spectrum list
          metadata=rbind.fill(metadata,metadata1) # Update metadata
          NN=NN+LL1
        }}
      }
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

  if (write_files){
    writeMGF2(library,output_library)
    write.table(library$metadata,paste0(output_library,".txt"),col.names = T,row.names=F,dec=".",sep="\t")}
  return(list(complete = library, current = library_current))
}
