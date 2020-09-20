#' Writing the spectral library to an open format
#'
#' The function writes the output of library_generator into mgf, msp or rdata file.
#'
#' @param output_library List that contains at list 2 elements.
#' \itemize{
#'   \item{sp:}{ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity}
#'   \item{metadata:}{ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on actually-detected scans. Following five columns are added: FILENAME, MSLEVEL, TIC, MASS_DEV, SCANNUMBER and SCANS}
#' }
#' @param con Character. Name of the output library, the file extension must be mgf, msp or rdata
#' @param type Character. "Complete" if the entire library is exported. "Consensus" if the consensus library is exported.
#'
#' @examples
#'
#' data(DRUG_THERMO_LIBRARY)
#' library_writer(library2,"library_V2_bis.mgf")
#'
#' @importFrom tools file_ext
#' @importFrom stringr str_remove
#'
#' @export
#' 
library_writer<-function(output_library, con = "output_library.mgf", type = "complete"){

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  ########################
  ### Check output format#
  ########################
  
  output_format = tolower(file_ext(con))
  
  if (!output_format %in% c("msp", "rdata", "mgf")){
      stop("The output library must be mgf, msp or RData format!")
  }
  
  type = type[1]
  if (!type %in% c("complete", "consensus")){
    stop("Please choose the exported library type between complete and consensus!")
  }
  
  if (is.null(output_library$consensus)){
    type = "complete"
    message("No consensus library is available! The entire library is exported!")
  }
  
  ################
  ###Read output##
  ################
  
  output_library = library_reader(output_library)
  
  if (type == "complete"){output_library = output_library$complete}
  if (type == "consensus"){output_library = output_library$consensus}
  
  #############
  ### Output ##
  #############
  
  if (output_format=="mgf"){writeMGF2(output_library, con = con)}
  if (output_format=="rdata"){save(output_library, file = con)}
  if (output_format=="msp"){writeMSP2(output_library, con = con)}
}

######################
### Internal function:
######################

writeMGF2 <- function(input_library, con) {
  
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  con <- file(description = con, open = "wt")
  on.exit(close(con))
  
  metadata = input_library$metadata
  splist = input_library$sp
  N=nrow(metadata)
  C=ncol(metadata)
  labels=colnames(metadata)
  for (i in 1:N) {
    .cat("\nBEGIN IONS\n")
    for (j in 1:C){
      .cat(labels[j],"=",as.character(metadata[i,j]),"\n")}
    sp=splist[[i]]
    .cat(paste(sp[,1],"\t",sp[,2], collapse = "\n"))
    .cat("\nEND IONS\n")
  }
}

writeMSP2 <- function(library, con) {
  
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  #######################
  # Conversion metadata:#
  #######################
  
  metadata = library$metadata
  splist = library$sp
  N1 = nrow(metadata)
  
  new_metadata = metadata
  
  ind = which(colnames(new_metadata)=="ID")
  colnames(new_metadata)[ind] = "Name"
  
  if ("SMILES" %in% colnames(new_metadata)){
    new_metadata$InChIKey = sapply(new_metadata$SMILES, inchikey_generator)
  } else {new_metadata$InChIKey = "N/A"}
  
  ind = which(colnames(new_metadata)=="ADDUCT")
  colnames(new_metadata)[ind] = "Precursor_type"
  
  ind = which(colnames(new_metadata)=="MSLEVEL")
  colnames(new_metadata)[ind] = "Spectrum_type"
  
  ind = which(colnames(new_metadata)=="PEPMASS")
  colnames(new_metadata)[ind] = "PrecursorMZ"
  
  ind = which(colnames(new_metadata)=="IONMODE")
  colnames(new_metadata)[ind] = "Ion_mode"
  temp = new_metadata[,ind]
  temp[temp=="Positive"] = "P"
  temp[temp=="Negative"] = "N"
  new_metadata[,ind] = temp
  
  ##########
  # Output:#
  ##########
  
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  con <- file(description = con, open = "wt")
  on.exit(close(con))
  
  for (i in 1:N1){
    
    ind = which(colnames(new_metadata)=="Name")
    .cat("Name: ",as.character(new_metadata[i,ind]),"\n")
    
    ind = which(colnames(new_metadata)=="InChIKey")
    .cat("InChIKey: ",as.character(new_metadata[i,ind]),"\n")
    
    ind = which(colnames(new_metadata)=="Precursor_type")
    .cat("Precursor_type: ",as.character(new_metadata[i,ind]),"\n")
    
    ind = which(colnames(new_metadata)=="Spectrum_type")
    .cat("Spectrum_type: ",as.character(new_metadata[i,ind]),"\n")
    
    ind = which(colnames(new_metadata)=="PrecursorMZ")
    .cat("PrecursorMZ: ",as.character(new_metadata[i,ind]),"\n")
    
    ind = which(colnames(new_metadata)=="Ion_mode")
    .cat("Ion_mode: ",as.character(new_metadata[i,ind]),"\n")
    
    .cat("Comments: \n")
    
    sp=splist[[i]]
    .cat("Num Peaks: ",nrow(sp),"\n")
    .cat(paste(sp[,1]," ",sp[,2], collapse = "\n"))
    
    .cat("\n")
    .cat("\n")
    
  }
}

#########
# Extra:#
#########

inchikey_generator<-function(smiles){
  
  options(stringsAsFactors = F)
  options(warn = -1)
  InChIKey = "N/A"
  
  base_url = "https://cactus.nci.nih.gov/chemical/structure/"
  url1 = paste0(base_url, smiles, "/InchiKey")
  temp1 = try(as.character(readLines(url1)[1]), silent = T)
  if (class(temp1)!="try-error"){InChIKey = str_remove(temp1, "InChIKey=")}
  
  return(InChIKey)
}
