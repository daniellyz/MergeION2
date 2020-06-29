#' Writing metadata and a list of spectra to mgf files
#'
#' The function writes a dataframe of metadata and a list of spectra to a msp file
#'
#' @param library List that contains 2 elements.
#' \itemize{
#'   \item{sp:}{ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity}
#'   \item{metadata:}{ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on actually-detected scans. Following five columns are added: FILENAME, MSLEVEL, TIC, MASS_DEV, SCANNUMBER and SCANS}
#' }
#' @param con Name of the output library, the file extension must be msp
#'
#' @examples
#'
#' data(DRUG_THERMO_LIBRARY)
#'
#' # Add new metadata "RESOLUTION = HIGH" to all scans:
#' library2$metadata$RESOLUTION = "HIGH"
#' # Write into a new mgf file:
#' writeMSP(library2,"library.msp")
#'
#' @importFrom tools file_ext
#' @importFrom stringr str_remove
#' 
#' @export
#'
writeMSP2 <- function(library, con) {

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  if (is.list(library)){
    if (length(library)==2 & "complete" %in% names(library)){
      library = library$complete
    }
  if (length(library)!=2 || (!is.list(library$sp)) || !is.data.frame(library$metadata)){
    stop("Please make sure your input library is a valid output of library_generator()!")
  }}

  if (is.character(con)){
    if (file_ext(con)!="msp"){
      stop("The file extension of your input library must be msp!")
    }} else {
      stop("The input must be the name of the msp file!")}

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

############
# Internal:#
############

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