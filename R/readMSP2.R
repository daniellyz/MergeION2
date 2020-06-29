#' Extracting metadata and a list of spectra from msp file
#'
#' The function extracts a dataframe of metadata and a list of spectra from msp file
#'
#' @param con Name of the input library, the file extension must be msp
#'
#' @return
#' \itemize{
#'   \item{sp:}{ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity}
#'   \item{metadata:}{ Data frame containing metadata of extracted scans.}
#' }
#' 
#' @importFrom tools file_ext
#' @importFrom plyr rbind.fill
#' 
#' @export
#'
readMSP2<-function(con){
  
  options(stringsAsFactors = FALSE)
  options(warn=-1)
  
  if (is.character(con)){
    if (file_ext(con)!="msp"){
      stop("The file extension of your input library must be msp!")
    }} else {
      stop("The input must be the name of the msp file!")}
  
  #####################################
  ### Reading from spectral library:###
  #####################################
  
  msp <- file(con,open="r")
  msp_lines <- readLines(msp,warn = F)
  msp_lines = msp_lines[msp_lines!=""]
  
  begin_ind = which(grepl("Name:",msp_lines))
  #end_ind = which(grepl("Num Peaks:",msp_lines))
  
  N = length(begin_ind)
  begin_ind = c(begin_ind, length(msp_lines)+1)

  metadata = c() # New metadata
  splist = list()

  for (i in 1:N){
    
    # Find metadata ranges:
    checked = FALSE
    ki = begin_ind[i]
    while (!checked){
      ki = ki + 1
      checked = !is.na(as.numeric(strsplit(msp_lines[ki], " ")[[1]][1]))
    }
    end_ind = ki - 1
  
    metadata_ind = begin_ind[i]:end_ind
    metadata_range = msp_lines[metadata_ind]
    NI = length(metadata_ind)
  
    new_metadata_names = sapply(metadata_range, function(x) strsplit(x,": ")[[1]][1])
    new_metadata_names = as.character(new_metadata_names)
    new_metadata = sapply(metadata_range, function(x) strsplit(x,": ")[[1]][2])
    new_metadata = as.data.frame(t(new_metadata))
    colnames(new_metadata) = new_metadata_names
    
    metadata = rbind.fill(metadata, new_metadata)
  
    spectrum_ind = (end_ind+1):(begin_ind[i+1]-1)
    nb_peak = length(spectrum_ind)
    spectrum_range = msp_lines[spectrum_ind]
    mzlist = sapply(spectrum_range,function(x) strsplit(x," ")[[1]][1])
    intlist = sapply(spectrum_range,function(x) strsplit(x," ")[[1]][2])
    spectrum = cbind(as.numeric(mzlist),as.numeric(intlist))
    splist[[i]] = spectrum
  }

  ####################
  ### Return results##:
  ####################

  rownames(metadata) = NULL
  metadata= data.frame(metadata)
  
  ind = which(colnames(metadata)=="Name")
  if (length(ind)==1){colnames(metadata)[ind] = "ID"}
  
  ind = which(colnames(metadata)=="Precursor_type")
  if (length(ind)==1){colnames(metadata)[ind] = "ADDUCT"}
  
  ind = which(colnames(metadata)=="Spectrum_type")
  if (length(ind)==1){colnames(metadata)[ind] = "MSLEVEL"}
  temp = metadata[,ind]
  temp[temp=="MS2"] = 2
  temp[temp=="MS1"] = 1
  metadata[,ind] = temp
  
  ind = which(colnames(metadata)=="PrecursorMZ")
  if (length(ind)==1){colnames(metadata)[ind] = "PEPMASS"}
  
  ind = which(colnames(metadata)=="Ion_mode")
  if (length(ind)==1){colnames(metadata)[ind] = "IONMODE"}
  temp = metadata[,ind]
  temp[temp=="P"] = "Positive"
  temp[temp=="N"] = "Negative"
  metadata[,ind] = temp
  
  ind = which(colnames(metadata)=="Comments")
  if (length(ind)==1){metadata = metadata[,-ind]}
  
  return(list(metadata = metadata, sp = splist))
}
