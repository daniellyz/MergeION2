#' Extracting metadata and a list of spectra from mgf file
#'
#' The function extracts a dataframe of metadata and a list of spectra from mgf file or
#'
#' @param con Name of the input library, the file extension must be mgf
#'
#' @return
#' \itemize{
#'   \item{sp:}{ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity}
#'   \item{metadata:}{ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on actually-detected scans. Following five columns are added: FILENAME, MSLEVEL, TIC, MASS_DEV, SCANNUMBER and SCANS}
#' }
#'
#' @importFrom tools file_ext
#' @importFrom MSnbase fData readMgfData
#'
#' @export
#'
readMGF2<-function(con){


  options(stringsAsFactors = FALSE)
  options(warn=-1)

  if (is.character(con)){
    if (file_ext(con)!="mgf"){
      stop("The file extension of your input library must be mgf!")
    }} else {
      stop("The input must be the name of the mgf file!")}

  #####################################
  ### Reading from spectral library:
  #####################################

  db = readMgfData(con, verbose = FALSE)
  metadata = fData(db)
  N = nrow(metadata)

  # From a MSnBase object to a list of spectra m/z intensity
  spectrum_list=list()
  for (i in 1:N){
    spectrum_list[[i]]=cbind(db[[i]]@mz,db[[i]]@intensity)
  }

  ####################
  ### Return results:
  ####################

  library = list()
  if (N>0){
    library$sp = spectrum_list
    library$metadata = metadata
  }
  return(library)
}
